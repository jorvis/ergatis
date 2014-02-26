#!/usr/bin/env perl

=head1 NAME

update_ec_numbers.pl - Will read in a current enzyme.dat file and update
    annotated EC numbers with suggestions from the file. Does not
    change partial EC numbers.

=head1 SYNOPSIS

 USAGE: update_ec_numbers.pl
       --input_ec_dat=/path/to/some/enzyme.dat
       --username=dbuser
       --password=dbpass
       --database_file=/path/to/db.txt
       --server=dbserver
       --no_change=1
     [ --password_file=/path/to/password.txt
       --log=/path/to/file.log
       --debug=3
       --help
     ]

=head1 OPTIONS

B<--input_ec_dat,-i>
    enzyme.dat file. (ftp://ftp.expasy.org/databases/enzyme/enzyme.dat)

B<--username,-u>
    Database username

B<--password,-p>
    Database password.  One of --password or --password_file is required
    
B<--password_file,-P>
    Password stored in a text file.  Useful for keeping password non-visible.

B<--database_file,-d>
    Either a comma-separated or tab-separated file containing the following columns:
    1) Name of Database
    2) Locus_tag ID prefix
    3) Path to a curate_common_names rules file
    4+) Any other DB related information (to come later)
    
B<--server,-s>
    Database hostname

B<--no_change,-n>
    Toggles between 0 (will change) and 1 (will not change). Will print out report but not actually change anything 
    in the database

B<--log,-l>
    Logfile.

B<--debug,-D>
    1,2 or 3. Higher values more verbose.

B<--help,-h>
    Print this message

=head1  DESCRIPTION

 Looks at enzyme.dat file (ftp://ftp.expasy.org/databases/enzyme/enzyme.dat)
 and makes sure all the EC numbers are current. Will delete whole EC numbers it 
 cannot find and will change EC numbers which it finds under Formerly EC tags
 in enzyme.dat file (basically, updating obsolete EC numbers)
 
=head1  INPUT
    Enzyme.dat (ftp://ftp.expasy.org/databases/enzyme/enzyme.dat)

=head1 OUTPUT
    None. Logfile. Report to STDOUT

=head1  CONTACT

    Kevin Galens
    kgalens@gmail.com

=cut

use strict;
use warnings;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Pod::Usage;
use DBI;
use Data::Dumper;

############# GLOBALS AND CONSTANTS ################
my $debug = 1;
my ($ERROR, $WARN, $DEBUG) = (1,2,3);
my $logfh;
my $dbh;
my $no_change = 0;
my %saved_queries;
my %report = ("removed"=>0, "changed"=>0, "genes"=>{});
my $password;
my $ec_rules;
my @undesired_ec = ();
####################################################

my %options;
my $results = GetOptions (\%options,
                          "username|u=s",
                          "password|p=s",
                          "password_file|P=s",
                          "database_file|d=s",
                          "server|s=s",
                          "input_ec_dat|i=s",
                          "no_change|n=i",
                          "log|l=s",
                          "debug|D=s",
                          "help|h"
                           );

&check_options(\%options);

if (defined $ec_rules) {
   &_log($DEBUG, "Parsing EC rules file $ec_rules");
   open RULES, $ec_rules or die ("Cannot run ec_rules file $ec_rules for reading: $!\n");
   @undesired_ec = <RULES>;
   close RULES;
}

&_log($DEBUG, "Parsing $options{'input_ec_dat'}" );
my $enzymes = &parse_dat_file( $options{'input_ec_dat'} );

&_log($DEBUG, "Getting transcripts..");
my $transcripts = &get_transcripts();

GENE:
    foreach my $t ( @{$transcripts} ) {
        my ($tid, $name) = @{$t};
        my $ecs = &get_ec_numbers( $tid );

      EC_NUMBER:
        foreach my $e ( @{$ecs} ) {
            # sadkins code to remove  particular obsolete ec nums from database
            my $found = 0;
            foreach my $obs (@undesired_ec) {
                chomp $obs;	# Did not chomp when parsing rules file
                if ($e->[0] eq $obs) {	# if our database ec matches a bad ec
                  $found = 1;
                  &_log($WARN, "Remove from $name: $e->[0] as the EC number is in rules file");
                  &remove_ec_number( $e->[1], $e->[0], $name );	# remove entry from db
                  last;
                }
            }
            next if ($found);	# Do not want to have any more operations on this EC number so go to next
            # end sadkins code

            ## skip if it isn't a whole ec number
            next EC_NUMBER if( $e->[0] =~ /-/ );
            if( exists( $enzymes->{'obsolete'}->{$e->[0]} ) ) {
                ## is it obsolete (i.e. can we change it?)
                &_log($WARN, "Change $name: $e->[0] to $enzymes->{'obsolete'}->{$e->[0]}");
                &change_ec_number( $enzymes->{'obsolete'}->{$e->[0]}, $e->[1], $name );
            } elsif( !exists( $enzymes->{'current'}->{$e->[0]} ) ) {
                &_log($WARN, "Remove from $name: $e->[0]");
                &remove_ec_number( $e->[1], $e->[0], $name );
            } 

        }
    }

print "No changes were actually made to database\n" if( $no_change );
print "Removed: ".$report{'removed'}."\n";
print "Updated: ".$report{'changed'}."\n";
print "Genes Affected: ".scalar(keys %{$report{'genes'}} )."\n";
$dbh->commit();
$dbh->disconnect();
################################### SUBS ###################################
sub remove_ec_number {
    my ($fcid, $ec, $name) = @_;

    $report{'removed'}++;
    $report{'genes'}->{$name} = 1;

    unless( $options{'no_change'} ) {
        my $query = 
            "DELETE FROM feature_cvterm WHERE feature_cvterm_id = ?";
        my $sth = $dbh->prepare( $query );
        $sth->execute( $fcid );
        $dbh->commit();
    } 
}

sub change_ec_number {
    my ($new, $fcid, $name) = @_;

    $report{'changed'}++;
    $report{'genes'}->{$name} = 1;

    unless( $options{'no_change'} ) {
        my $cvterm_id = &get_cvterm_id_by_ec( $new );
        my $query = 
            "UPDATE feature_cvterm SET cvterm_id = ? ".
            "WHERE feature_cvterm_id = ?";
        my $sth = $dbh->prepare( $query );
        $sth->execute( $cvterm_id, $fcid );
        $dbh->commit();
    }
}

sub get_cvterm_id_by_ec {
    my ($ec) = @_;
    my $cvterm_id;
    my $query = 
        "SELECT c.cvterm_id FROM cvterm c, cv, dbxref d, cvterm_dbxref cd ".
        "WHERE cd.cvterm_id = c.cvterm_id AND cd.dbxref_id = d.dbxref_id ".
        "AND cv.cv_id = c.cv_id AND cv.name = 'EC' ".
        "AND d.accession = '$ec'";

    my $r = &do_query($query);
    &_log($ERROR, "Expected 1 result from query, returned ".scalar(@{$r})) if( @{$r} > 1 );
    
    ## this means that the ec number doesn't exist
    if( @{$r} == 0 ) {
        $cvterm_id = &insert_ec_number( $ec );
    } else {
        $cvterm_id = $r->[0]->[0];
    }

    return $cvterm_id;
}

sub get_next_id {
    my ($table, $col) = @_;
    my $query = 
        "SELECT max( $col ) FROM $table";
    my $r = &do_query( $query );
    my $num = $r->[0]->[0];
    ++$num;
}

sub insert_ec_number {
    my ($ec) = @_;
    
    &_log($ERROR, "Could not find [$ec] in enzyme lookup when inserting")
        unless( exists( $enzymes->{'current'}->{$ec} ) );
    my $ec_info = $enzymes->{'current'}->{$ec};

    ## grab the db_id and cv_id
    my $db_id = &get_id_by_name("EC", "db");
    my $cv_id = &get_id_by_name("EC", "cv");

    ## get next cvterm_id, dbxref_id and cvterm_dbxref_id
    my $cvterm_id = &get_next_id( 'cvterm', 'cvterm_id' );
    my $dbxref_id = &get_next_id( 'dbxref', 'dbxref_id' );
    my $cd_id = &get_next_id("cvterm_dbxref", "cvterm_dbxref_id");

    ## insert cvterm
    &do_insert_query( "INSERT INTO cvterm (cvterm_id, cv_id, name, definition, is_obsolete, is_relationshiptype) ".
                      "VALUES( ?, ?, ?, ?, 0, 0 )", $cvterm_id, $cv_id, $ec_info->{'name'}, $ec_info->{'def'} );

    ## insert dbxref
    &do_insert_query( "INSERT INTO dbxref (dbxref_id, db_id, accession, version) ".
                      "VALUES (?, ?, ?, 'current')", $dbxref_id, $db_id, $ec);

    ## insert cvterm_dbxref
    &do_insert_query( "INSERT INTO cvterm_dbxref (cvterm_dbxref_id, cvterm_id, dbxref_id, is_for_definition) ".
                      "VALUES (?, ?, ?, 0)", $cd_id, $cvterm_id, $dbxref_id );
    

    ## commit changes
    $dbh->commit();

    return $cvterm_id;
}

sub do_insert_query {
    my ($query, @args) = @_;
    my $sth = $dbh->prepare( $query );
    $sth->execute( @args );
}

sub get_id_by_name {
    my ($name, $table) = @_;
    my $query = 
        "SELECT ${table}_id FROM $table WHERE name = ?";
    my $r = &do_query( $query, $name );
    $r->[0]->[0];
}

sub get_ec_numbers {
    my ($tid) = @_;

    my $query = 
        "SELECT d.accession, fc.feature_cvterm_id ".
        "FROM feature f, feature_cvterm fc, cvterm c, dbxref d, cv, cvterm_dbxref cd ".
        "WHERE f.feature_id = fc.feature_id ".
        "AND fc.cvterm_id = c.cvterm_id ".
        "AND c.cv_id = cv.cv_id ".
        "AND cv.name = 'EC' ".
        "AND c.cvterm_id = cd.cvterm_id ".
        "AND cd.dbxref_id = d.dbxref_id ".
        "AND f.feature_id = ?";
    &do_query( $query, $tid );
}

sub get_transcripts {
    my $query =
        "SELECT t.feature_id, t.uniquename ".
        "FROM feature t, cvterm c ".
        "WHERE t.type_id = c.cvterm_id ".
        "AND c.name = 'transcript'";
    my $results = &do_query( $query );
    return $results;
}

sub parse_dat_file {
    my ($file) = @_;
    open(IN, "< $file") or die("Could not open $file: $!");

    my %data = ('current'=>{},'obsolete'=>{});
    my ($current, $transfer, $def, $name, $deline);
    while( my $line = <IN> ) {
        chomp($line);
        if( $line =~ m|^//| ) {
            
            if( defined( $current )  ) {
                if( $deline =~ /Transferred entry:\s+(.*)/ ) {
                    #Sometimes we have multiple entries here:
                    #example: "1.5.3.13, 1.5.3.14, 1.5.3.15, 1.5.3.16 and 1.5.3.17"
                    my $st = $1;
                    my @r = split(/,\s/, $st);
                    my @ecs;
                    map { 
                        my @a = split(/\sand\s/, $_); 
                        map { $_ =~ s/\.$// } @a;
                        push(@ecs, @a); 
                    } @r;

                    my $tran;
                    if( @ecs > 1 ) {
                        $tran = &find_common_parent( @ecs );
                    } else {
                        $tran = $ecs[0]
                    }
                    $transfer = $tran;
                    $transfer =~ s/\.$//;
                } else {
                    $name = $deline." [$current]";
                }
                
                ## we should store current
                if( $transfer ) {
                    $data{'obsolete'}->{$current} = $transfer;
                } else {
                    $data{'current'}->{$current} = {
                        'name' => $name,
                        'def' => $def
                        };
                }
            }

            ## reset the variables
            map { undef($_); } ($current, $transfer, $def, $name, $deline);

        } elsif( $line =~ /^ID\s+(.*)/ ) {
            $current = $1;
        } elsif( $line =~ /^DE\s+(.*)/ ) {
            $deline = " ".$1;
            $deline =~ s/^\s+//;
        } elsif( $line =~ /^CA\s+(.*)/ ) { 
            $def = $1;
        }
    }
    close(IN);
    return \%data;
}

sub find_common_parent {
    my (@nums) = @_;
    
    my %ec_numbers;
    map { 
        my @t = split(/\./, $_);
        $ec_numbers{$_} = \@t;
    } @nums;
    
    my @arr = qw(- - - -);
    POS:    
    for(my $i = 0; $i < 4; $i++ ) {
        my $cur_value;
        EC_NUMBER:
        foreach my $ecnum ( keys %ec_numbers ) {
            my $val = $ec_numbers{$ecnum}->[$i];
            if( !defined( $cur_value ) ) {
                $cur_value = $val;
            } else {
                if( $val =~ /and/ || $cur_value =~ /and/ ) {
                    print Dumper( \@nums );
                    print "Val: $val, Cur value $cur_value\n";
                    die("Shouldn't happen!");
                }
                if( $val != $cur_value ) {
                    last POS;
                }
            }
        }
        $arr[$i] = $cur_value;
    }
    my $retval = join(".", @arr);
    return $retval;
}

sub do_query {
    my ($query, @args) = @_;
    &_log($ERROR, "No query defined") unless( defined( $query ) );
    my $sth;
    if( exists( $saved_queries{$query} ) ) {
        $sth = $saved_queries{$query};
    } else {
        $sth = $dbh->prepare( $query );
        $saved_queries{$query} = $sth;
    }
    $sth->execute( @args );
    my $r = $sth->fetchall_arrayref();
    #die("Query returned 0 results: $query [@args]") unless( @{$r} );
    return $r;
}


############################# UTIL ###############################
sub check_options {

   my $opts = shift;

   if( $opts->{'help'} ) {
       &_pod;
   }

   if( $opts->{'log'} ) {
       open( $logfh, "> $opts->{'log'}") or die("Can't open log file ($!)");
   }

   foreach my $req ( qw(input_ec_dat username database_file server) ) {
       &_log($ERROR, "Option $req is required") unless( $opts->{$req} );
   }

   open DB, $options{'database_file'} or die "Cannot open database_file for reading: $!\n";
   my $line = <DB>;
   chomp $line;
   my @metadata = split(/\t/, $line);   
   my $database = $metadata[0];
   $ec_rules = $metadata[18] if (defined $metadata[18]);
   close DB;

    #Assign password to be read from the file if it exists.
    if (defined ($options{'password_file'} && -s $options{'password_file'} ) ) {
	  open PASS, $options{'password_file'} or die ("Cannot open password file ". $options{'password_file'} ." : $!\n");
	  print STDERR ("Password from file will take priority over --password option\n") if (defined $options{'password'});
	  $password= <PASS>;
	  chomp $password;
	  close PASS;
    } elsif (defined ($options{'password'}) ){
	  $password = $options{'password'};
    } else {
	  die("Neither a password or a password file were supplied.  Please supply one or the other");
    }

   $debug = $opts->{'debug'} if( $opts->{'debug'} );

   $dbh = &_connect( $opts->{'username'}, $password, $database,
                     $opts->{'server'} );

    if ( defined($opts->{'no_change'}) && ($opts->{'no_change'} != 0 && $opts->{'no_change'} != 1) ) {
   		die("The --no_change option must either be a 0 (will change) or a 1 (will not change).");
    }
   $no_change = 1 if( $opts->{'no_change'} == 1 );
}

sub _log {
   my ($level, $msg) = @_;
   if( $level <= $debug ) {
      print STDOUT "$msg\n";
   }
   print $logfh "$msg\n" if( defined( $logfh ) );
   exit(1) if( $level == $ERROR );
}


sub _connect {
  my ($user, $password, $db, $server) = @_;
  my $dbh;
  eval {
  $dbh = DBI->connect("DBI:mysql:$db:$server", "$user", "$password",
                       {
                                'RaiseError' => 1,
                                'AutoCommit' => 0,
                          } );
    };
    if( $@ ) {
        die("Could not connect to database ".DBI->errstr);
    }
    $dbh->do("use $db");
    return $dbh;
}

sub _pod {
    pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDERR} );
}
