#!/usr/local/bin/perl -w

=head1 NAME

wublast_btab_loader.pl - Loads wu-blast hits from btab into a legacydb

=head1 SYNOPSIS

    USAGE: wublast_btab_loader.pl OPTION LIST HERE

=head1 OPTIONS

B<--ev_type,-e>
    - The value to store in evidence..ev_type

B<--input_file,-f>
    - When giving a single file as input for loading.

B<--input_list,-i>
    - When giving a list of files as input for loading.

B<--search_db, -s>
    - Which db was searched against?

B<--project_db, -p>
    - name of project database

B<--user, -U>
    - user to log into the db server as

B<--password, -P>
    - password associated with the given userid

B<--log,-l>
    - path to log file

B<--delete_old,-x>
    - remove old evidences for ev_type on the encountered asmbl_ids

B<--min_match_length,-m>
    - filter hits on the match length prior to loading

B<--min_per_id,-n>
    - filter hits on the percent identity prior to loading

B<--min_per_cov,-o>
    - filter hits on the percent coverage prior to loading

B<--help,-h>
    - display this help text

B<--debug,-d>
    - log individual queries.

=head1  DESCRIPTION

Loads blast btab files into a legacydb. 

=head1  INPUT

Takes in a btab file or a list of such files.  Filenames must be in the form:
/path/[db].assembly.[asmbl_id].wu-blast[n|p].btab

=head1  OUTPUT

Output consists of a log, although if successful, a legacy db will contain the 
results once the script is finished running.

=head1  CONTACT

    Jason Inman
    jinman@tigr.org

=begin COMMENT
## legal values for status are active, inactive, hidden, unstable
    status: active
    keywords: keywords to search for your script here
=end COMMENT

=cut

use strict;
use warnings;
$|++;

use Pod::Usage;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use DBI;

my @handles;
my @inputFiles;
my $project_db;
my $search_db;
my $db_id;
my $ev_type;
my $min_match_length = 0;
my $min_per_id = 0;
my $min_per_cov = 0;
my $delete_old = 0;
my $user = 'egc';
my $password = 'egcpwd';

my %options = ();
my $debug = 0;
my $results = GetOptions(\%options,
                        'input_file|f=s',
                        'input_list|i=s',
                        'ev_type|e=s',
                        'search_db|s=s',
                        'project_db|p=s',
                        'log|l=s',
                        'delete_old|x',
                        'min_match_length|m=f',
                        'min_per_id|n=f',
                        'min_per_cov|o=f',
                        'debug=i',
                        'help|h',
                        'user|U=s',
                        'password|P=s',
                        ) || &_pod;

&check_parameters(\%options);

my $dbh = DBI->connect("dbi:Sybase:server=SYBTIGR","egc","egcpwd",{PrintError => 0});
&_die ("Can't connect to SYBTIGR as egc") unless defined $dbh;
$dbh->do("use $project_db") || &_die ("Can't use $project_db: $!");
my $sth;

my $BTAB_EVIDENCE_INSERTION = "INSERT evidence (feat_name, m_lend, m_rend, end5, end3, change_log,
                                      save_history, accession, assignby, curated, ev_type,
                                      method, per_id, per_sim, score, pvalue, rel_end5,
                                      rel_end3, date, db)
                           VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?,
                                   ?, ?, ?, ?, 0, 0, getdate(), ?)";


my $insert_btab_evidence = $dbh->prepare($BTAB_EVIDENCE_INSERTION);

# Get this ready if we're going to be deleting:
my $del_handle if $delete_old;
my $find_old   if $delete_old;
if ($delete_old) {

    my $FIND_OLD_EVS = "SELECT e.id FROM evidence e, common..search_dbs s
                        WHERE e.ev_type = ? AND convert(numeric(10),e.db) = s.id
                        AND s.name = ? and e.feat_name like ?";

    $find_old   = $dbh->prepare($FIND_OLD_EVS);

    $del_handle = $dbh->prepare("DELETE FROM evidence WHERE id = ?");    

}

## Ensure the search db is in common..search_dbs
&check_search_db;

# Go through each input file and handle it
foreach my $file (@inputFiles) {

    &_log ("Working on $file\n");

    my $asmbl_id = &get_asmbl_id($file);

    # Handle deletions if necessary
    if ($delete_old) {
        &delete_existing_evidence($asmbl_id);
    }

    # load. 
    &load_btab($file,$ev_type, $search_db, $asmbl_id, $min_match_length,
                                                    $min_per_id, $min_per_cov);

}

close $handles[0];

exit(0);

######### Subprocedures
sub check_parameters {
# get all our options in order.
    
    my $options = shift;

    my $error = '';

    if ($options{'help'}) {
        &_pod;
        exit;
    }

    # Set up the input file list
    if ($options{'input_list'}) {
        if (-e $options{'input_list'}) {
            open (IN, "$options{'input_list'}") ||
                die ("Couldn't open $options{'input_list'}");
            @inputFiles = <IN>;
            close(IN);
            chomp foreach @inputFiles;
        } else {
            $error .= "$options{'input_list'} doesn't exist\n";
        }
    } elsif ($options{'input_file'}) {
        if (-e $options{'input_file'}) {
            push @inputFiles, $options{'input_file'};
        } else {
            $error .= "$options{'input_file'} doesn't exist\n";
        }
    } else {
        $error .= "Must use either --input_list or --input_file\n";
    }

    if ($options{'project_db'}) {
        $project_db = $options{'project_db'};
    } else {
        $error .= "project_db must be specified.\n";
    }

    if ($options{'search_db'}) {
        $search_db = $options{'search_db'};
    } else {
        $error .= "search_db must be specified\n";
    }

    if ($options{'log'}) {
        open (LOGFILE, "> $options{'log'}") ||
            die ("Can't write to $options{'log'}: $!");
        push @handles, *LOGFILE;
    } else {
        push @handles, *STDOUT;
    }

    if ($options{'ev_type'}) {
        $ev_type = $options{'ev_type'};
    } else {
        $error .= "must specify ev_type\n";
    }

    if ($options{'min_match_length'}) {
        $min_match_length = $options{'min_match_length'};
    }

    if ($options{'min_per_id'}) {
        $min_per_id = $options{'min_per_id'};
    }

    if ($options{'min_per_cov'}) {
        $min_per_cov = $options{'min_per_cov'};
    }

    if ($options{'debug'}) {
        $debug = $options{'debug'};
    }

    if ($options{'delete_old'}) {
        $delete_old = $options{'delete_old'};
    }

    if ($options{'user'}) {
        $user = $options{'user'};
        if ($options{'password'}) {
            $password = $options{'password'};
        } else {
            print "Password for user $user: ";
            system('stty', '-echo');
            chomp($password = <STDIN>);
            system('stty', 'echo');
            print "\n";
        }
    }

    if ($error) {
        &_die("$error");
    }
    

}

sub check_search_db {
## Before parsing the btabs, make sure the search db has an id in common..search_dbs

    my $query = "SELECT id FROM common..search_dbs
                 WHERE name LIKE \"$search_db\"
                 AND iscurrent = 1";

    $sth = $dbh->prepare($query);
    $sth->execute();

    my @results = $sth->fetchrow();

    unless ($db_id = $results[0]) {
        my $warn = "The search_db $search_db is not listed as a db in common..search_dbs\n"
        .  "and therefore, will not be loaded into the annotation database: $project_db\n";
        &_die($warn);
    }
    &_log ("db id for $search_db is $db_id\n") if $debug;
}

sub delete_existing_evidence {
## Delete old search data for this assembly against this dba

    my $curr_asmbl_id = shift;    
    print "Deleting old search evidence for ev_type: $ev_type and database: $search_db " .
          "on assembly: $curr_asmbl_id\n";

    $find_old->execute($ev_type,$search_db,"$curr_asmbl_id\.%");
    
    my $results_ref = $find_old->fetchall_arrayref();
    my $num_found =  scalar(@$results_ref);

    if ($num_found == 0) {

        &_log("\tNo results found to delete\n");

    } else {

        &_log("\tRemoving $num_found evidence records\n");

        foreach my $row (@$results_ref) {

            my ($id) = @$row;
            &_log("\t\tRemoving evidence with id: $id\n") if $debug;
            ## DELETE STATEMENT GOES HERE:
            $del_handle->execute($id);

        }
    }

}

sub load_btab {
# Given a btab file:
# Load the btab.
    my ($btabfile, $ev_type, $search_db, $asmbl_id, $min_match_length, $min_per_id, $min_per_cov) = @_;

    &_log("\tLoading new search results\n");

    open (BTAB, $btabfile) || die "cant open the file $btabfile";
    my $line;
    my $line_count = 0;
    while ($line = <BTAB>) {

##IS THIS NEEDED?  I don't think so.  jinman##        $line =~ s/\tbac://; #bhaas hack.
        next unless ($line =~ /\w/);
        #print "$line";
        chomp ($line);
        my @btab = split(/\t/, $line);

#  0 : zma1.assembly.8844   #### query id
#  1 :                      #### empty in wu-blast btab
#  2 : 9396                 #### query sequence length
#  3 : BLASTN               #### method
#  4 : MF                   #### searchdatabase
#  5 : OG0EG31TV            #### accession
#  6 : 7359                 #### end5
#  7 : 8305                 #### end3
#  8 : 1                    #### match left
#  9 : 947                  #### match right
#  10 : 100.00              #### per_id
#  11 : 100.00              #### per_sim
#  12 : 4735                #### score
#  13 : 716.5               #### 
#  14 :                     #### 
#  15 :                     ####
#  16 : 1                   #### strand as numeral (-1 = minus, 1 = plus)
#  17 : Plus                #### strand as word
#  18 : 947                 #### match length
#  19 : 4.1e-207            #### p-value

##        my $method = $btab[3]; 
## I'm choosing to store the program name of this loader, as it seems to be what most 
## entries are using evidence.evtype for. However, I personally think something like "BLASTN"
## or the like would be more appropriate.
        my $method = $0;

        my $accession = $btab[5];
        my $end5 = $btab[6];
        my $end3 = $btab[7];
        my $m_lend = $btab[8];
        my $m_rend = $btab[9];
        my $per_id = $btab[10];
        my $per_sim = $btab[11];
        my $score = $btab[12];
        my $match_len = $btab[18];
        my $pvalue = $btab[19];
        unless ($pvalue) {
            $pvalue = "NULL";
        }

        next if ( ($min_per_id) && ($per_id < $min_per_id));

        my $per_cov = ( (abs($m_lend - $m_rend) + 1) / $match_len ) * 100;
        next if ($per_cov < $min_per_cov);

        my $query_match_length = abs ($end5 - $end3) + 1;
        my $db_match_length = abs ($m_lend - $m_rend) + 1;
        next unless ($db_match_length >= $min_match_length 
                        && $query_match_length >= $min_match_length);

&_log("Attempting load with these values:\n\t\"$asmbl_id.intergenic\", $m_lend, $m_rend,$end5, $end3, 0, 0, $accession, $user, 0, $ev_type, $method, $per_id, $per_sim,$score, $pvalue, 0, 0, getdate(), $db_id");

        $insert_btab_evidence->execute("$asmbl_id.intergenic", $m_lend, $m_rend, 
                                        $end5, $end3, 0, 0, $accession, $user, 0,
                                        $ev_type, $method, $per_id, $per_sim,
                                        $score, $pvalue, $db_id);
        
        $line_count++;

    }

    close (BTAB);
    if ($line_count) {
        &_log("Loaded $line_count records for $btabfile"); 
    } else {
       &_log("Nothing loaded for $btabfile"); 
    }
        
}


sub get_asmbl_id {
# Given a filename, return the portion corresponding to the asmbl_id

    my $filename = shift;

    if ($filename =~ /$project_db\.\w+\.(\d+).*btab$/) {
        return $1;
    } else {
        &_die ("Bad filename: Can't extract asmbl_id from: $filename\n");
    }

}

sub _pod {
#Print the help text
    pod2usage( {-exitval => 0, -verbose => 2 });
}

sub _log {
# send output to all filehandles (well, the one or two filehandles) in @handles
    my $msg = shift;
    chomp $msg;

    foreach (@handles) {
        print $_ "$msg\n";
    }

}

sub _die {
# Print a message to the log and then close it, then die.
    my $msg = shift;
    chomp $msg;

    &_log("$msg");
    close $handles[0];

    die "$msg\n";

}

sub do_sql {
# take in a dbi object, create a statement handler and return a reference to the 
# results

    my ($dbh, $query) = @_;

    my $sth = $dbh->prepare($query);
    $sth->execute();

    
}
