#!/usr/bin/env perl

=head1 NAME

remove_db_gene_syms_with_hypos.pl - Description

=head1 SYNOPSIS

 USAGE: remove_db_gene_syms_with_hypos.pl
       --database_file=/path/to/database.file
       --username=db_user1
       --password=db_password1
       --update=1
       [--password_file=/path/to/password/txt]
       --server=annot_db.someplace.org

=head1 OPTIONS

B<--database_file,-d>
    Either a comma-separated or tab-separated file containing the following columns:
    1) Name of Database
    2) Locus_tag ID prefix
    3) Path to a curate_common_names rules file
    4+) Any other DB related information (to come later)
    
B<--username,-u>
    The MySQL username of a user with SELECT permissions on --database.

B<--password,-p>
    The MySQL password for the username specified by --username on the MySQL server --server.

B<--password_file,-P>
	Password stored in a text file.  Useful for keeping password non-visible.

B<--server,-s>
    The MySQL server on which the database specified by the --database option resides.
    
B<--update,-a> 
	Default = 0. Set to non-zero to make changes to database. Will not change anything by default    
    
B<--log,-l>
    Optional.  Location of a log file to which error messages/warnings should be written.

B<--help,-h>
    Optional.  Print this documentation.

=head1  DESCRIPTION

This script removes, from the database, any entries with a gene symbol and a 'hypothetical protein' gene product name
Matches "hypothetical protein", "conserved hypothetical protein", "conserved domain protein", and "conserved protein"

=head1  CONTACT
    Shaun Adkins
    sadkins@som.umaryland.edu

=cut

use strict;
use warnings;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Pod::Usage;
use Data::Dumper;
use DBI;

############# GLOBALS AND CONSTANTS ################
my $debug = 1;
my ($ERROR, $WARN, $DEBUG) = (1,2,3);
my $logfh;

my $output_file;
my $DBH;
my %QUERIES = ();
my $password;
my $database;
my $dbh;
my $update = 0;
####################################################

my %options;
my $results = GetOptions (\%options,
                          "database_file|d=s",
                          "username|u=s",
                          "password|p=s",
                          "password_file|P=s",
                          "server|s=s",
                          "output_file|o=s",
                          "update|u=s",
                          "log|l=s",
                          "debug=s",
                          "help|h"
                          );

check_options(\%options);

$dbh = &_connect($options{'username'}, $password, $database, $options{'server'});
prepare_queries( \%QUERIES );  # prepare SQL commands to issue to the Chado database	
if (!$update) {	# Select instead of delete if we do not plan to update db
	$QUERIES{'select_gene_syms_with_hypos'}->execute();
	print STDOUT "There were ". scalar($QUERIES{'select_gene_syms_with hypos'}->rows) . " entries found in " . $database, "\n";
} else {
	$QUERIES{'remove_gene_syms_with hypos'}->execute();
	print STDOUT "There were ". scalar($QUERIES{'remove_gene_syms_with hypos'}->rows) . " entries removed from " . $database, "\n";
}
exit(0);

sub prepare_queries {
	# This is where the mySQL queries are prepared to be later executed in the scripts
	my ($queries) = @_;
	
	my $remove = "DELETE fp1 FROM featureprop as fp1 " . 
			"INNER JOIN featureprop as fp2 " .
			"WHERE fp1.feature_id = fp2.feature_id " .
			"AND fp1.type_id = 65 AND fp2.type_id = 68 " .
			"AND fp1.value IS NOT NULL " .
			"AND (fp2.value = 'hypothetical protein' " .
			"OR fp2.value = 'conserved hypothetical protein' " .
			"OR fp2.value = 'conserved domain protein' " .
			"OR fp2.value = 'conserved protein')";
			
	#cvterm_id 65 = gene_symbol
	#cvterm_id 68 = gene_product_name			
				
	$queries->{'remove_gene_syms_with hypos'} = $DBH->prepare($remove);
	
	my $select = "SELECT fp1 FROM featureprop as fp1 " . 
			"INNER JOIN featureprop as fp2 " .
			"WHERE fp1.feature_id = fp2.feature_id " .
			"AND fp1.type_id = 65 AND fp2.type_id = 68 " .
			"AND fp1.value IS NOT NULL " .
			"AND (fp2.value = 'hypothetical protein' " .
			"OR fp2.value = 'conserved hypothetical protein' " .
			"OR fp2.value = 'conserved domain protein' " .
			"OR fp2.value = 'conserved protein')";
	$queries->{'select_gene_syms_with_hypos'} = $DBH->prepare($remove);
}

sub check_options {
   my $opts = shift;
   if( $opts->{'help'} ) {
       &_pod;
   }

   if( $opts->{'log'} ) {
       open( $logfh, "> $opts->{'log'}") or die("Can't open log file ($!)");
   }

   $debug = $opts->{'debug'} if( $opts->{'debug'} );

   foreach my $req ( qw(database_file username server) ) {
       &_log($ERROR, "Option $req is required") unless( $opts->{$req} );
   }
   
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
   
    open DB, $options{'database_file'} or die ("Cannot open database file " . $options{'database_file'} . "$!\n");
   	my $line = <DB>;
   	chomp $line;
   	my $rest;
   	($database, $rest) = split(/,|\t/, $line, 2);
    close DB;
   
    $update = 1 if( $options{'update'} );
   
}

sub _connect {
  my ($user, $password, $db, $server) = @_;
  eval {
  $DBH = DBI->connect("DBI:mysql:$db:$server", "$user", "$password",
                       {
                                'RaiseError' => 1,
                                'AutoCommit' => 1,
                          } );
    };
    if( $@ ) {
        die("Could not connect to database ".DBI->errstr);
    }
    $DBH->do("use $db");
}

sub _log {
   my ($level, $msg) = @_;
   if( $level <= $debug ) {
      print STDOUT "$msg\n";
   }
   print $logfh "$msg\n" if( defined( $logfh ) );
   exit(1) if( $level == $ERROR );
}

sub _pod {
    pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDERR} );
}
