#!/usr/local/bin/perl

=head1 NAME
prokpipe_consistency_checks.pl
    Automated pipeline and database error detection
    Checks BSML files and DB for consistencies and either warns or errors out
	dependng on the failed check and its severity
    Used in the Prokaryotic Annotation Pipeline run by Ergatis

=head1 OPTIONS (more to be added later)
B<--input_list, -i>
    A list containing the paths of .bsml files

B<--output_dir, -o>
    Name of the output directory
    
B<--database, -d>
    Use a specified database instead of an input list
    
B<--host, -h>
    Host of the database if you are checking that
    
B<--user, -u>
    Username for the database
    
B<--pass, -p>
    Password for your username
    
B<--log,-l>
    Logfile.

B<--help>
    Print this message

=head1 SYNOPSIS

    For checking BSML lists:
    error_detector.pl -i path/to/files.bmsl.list -o path/to/output/dir
    
    For checking databases:
    prokpipe_consistency_checks.pl -d database -h host -u username -p password -o path/to/output/dir
    
=head1 CONTACT
    Shaun Adkins
    sadkins@som.umaryland.edu

=cut

use strict;
use warnings;

use File::Basename;
use FindBin qw($Bin);
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Pod::Usage;
use Prokpipe_Checks::BSMLChecker;
use Prokpipe_Checks::DBChecker;
use Ergatis::Logger;
use DBI;


##### Prototypes #####
sub check_parameters($);
sub _connect($$$$);
sub read_file($);
sub run_checks($);
sub run_checks_db();    #to be implemented later

my %options = ();
my $results = GetOptions (\%options,
		'output_dir|o=s',
		'input_list|i=s',
		'database|d=s',
                'user|u=s',
                'pass|p=s',
                'host|h=s',
                'log|l=s',
		'debug|b=s',
                'help') || pod2usage();

my $database_flag = 0;
my ($ERROR, $WARN, $DEBUG) = (1,2,3);
my $dbh;

my $database;
my $type;
my @paths;
my $val_file;

## Display documentation
if( $options{'help'} ){
    pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDERR} );
}

## Getting the log file
my $logfile = $options{'log'} || Ergatis::Logger::get_default_logfilename();
our $logger = new Ergatis::Logger('LOG_FILE'=>$logfile,
                'LOG_LEVEL'=>$options{'debug'});
$logger = $logger->get_logger();

check_parameters(\%options);

(!$database_flag) ? $val_file = $options{'output_dir'}."/db.validation.out" : $val_file = $options{'output_dir'}."/bsml.validation.out";
open VF, ">$val_file" or $logger->logdie("Could not open $val_file file for writing: $!\n");


(!$database_flag) ? run_checks(\@paths) : run_checks_db();

close VF;


##### Subroutines #####

# Subroutine to run various checks to make sure all is well
sub run_checks($) {
    my $paths = shift;
    my $state = 0;
    my $name = "";
    my $id = "";
    my $warn_flag = 0;	#currently have not decided what to do with this... may have separate one for warn and error
    my $error_flag = 0;	# Plan on killing component after running tests, to get all errors present
    my $length_flag = 0; # Flag is raised if the assembly length is less than the cutoff (1000 bp)
    
# Check for proper feature-type count in each file
# If a file contains none of a particular type, the module will return an error
    foreach my $file (@{$paths}){
	print VF ">$file\n";
	$state = $type->get_length($file);
	$length_flag = 1 if ($state);

	if (!$length_flag) {
        $state = $type->count_genes($file);
        $logger->logdie("ERROR, No genes present in $file")if ($state == $ERROR);

        $state = $type->count_cds($file);
        $logger->logdie("ERROR, No CDS domains present in $file")if ($state == $ERROR);

        $state = $type->count_polypeptides($file);
        $logger->logdie("ERROR, No polypeptides present in $file")if ($state == $ERROR);

        $state = $type->count_transcripts($file);
        $logger->logdie("ERROR, No transcripts present in $file")if ($state == $ERROR);

        $state = $type->count_exons($file);
        $logger->logdie("ERROR, No exons present in $file")if ($state == $ERROR);
	
	$type->count_tRNA($file);
	$warn_flag = 1 if ($type->divide_by_3($file));
#	$warn_flag = 1 if ($type->bad_gene_symbols($file));
	$warn_flag = 1 if ($type->duplicate_gene_symbols($file));
	$warn_flag = 1 if ($type->should_not_have_gs($file));
	$warn_flag = 1 if ($type->ec_check($file));
	$warn_flag = 1 if ($type->TIGR_role_check($file));
	$warn_flag = 1 if ($type->valid_start_stop_codons($file));
    	}

	print VF $type->return_counts();
	$type->_reset();
    }
    
    # Add more checks for the perl module to check as time passes
    
    print "Looks like everything is A-ok for loading!\n";

    print VF "The file $options{'input_list'} has no errors and is die for Chado loading\n";
    #This is more of a test output since every component has to have an output.
    #Once we know what checks we want to run, then the output will be adjusted accordingly
}

# Subroutine to run various checks on a db
sub run_checks_db() {
    print "TO BE IMPLEMENTED LATER\n";

    my $die = 0;
    
    $die = $type->no_GO_evidence;
    $logger->logdie("ERROR, GO Role Links found with no evidence") if ($die == $ERROR);
    $die = $type->no_GS_EC;
    $logger->logdie("ERROR, No gene symbol or EC number present in conserved hypotheical protein") if ($die == $ERROR);#change error later

}

# Subroutine to check if the supplied paramaters are correct
sub check_parameters($) {
    my $options = shift;

    my $req;

	## Confirm we are working with databases or not
    if (defined($options{'database'})) {
        $database_flag = 1;
        $database = $options{'database'};
        $type = "Prokpipe_Checks::DBChecker";    #Decide which perl module to use for the checks
        foreach $req ( qw(user pass database host output_dir) ) {
            $logger->logdie("ERROR: Option $req is required") unless( $options{$req} );
        }
        _connect( $options{'user'}, $options{'pass'}, $options{'database'}, $options{'host'} );
    } else {
        $type = "Prokpipe_Checks::BSMLChecker";
	## make sure the input list file exists and is readable
        foreach $req ( qw(input_list output_dir) ) {
            $logger->logdie("ERROR: Option $req is required") unless( $options{$req} );
        }
        if (! -e "$options{'input_list'}") {
        	$logger->logdie("The $options{'input_list'} input list file passed could not be read or does not exist");
        } else {
        	@paths = read_file($options{'input_list'});
        	foreach my $bsml (@paths) {
        	    chomp($bsml);
                next if ($bsml =~ /^\s*$/);	#ignore whitespace
    	    next if ($bsml =~ /^#/);	#ignore 
                if ($bsml =~ /\//g) {           #Make sure the list contains paths to files
                } else {
                    $logger->logdie("The $options{'input_list'} file is not a list file. Incorrect input");
                }
            }
        }
    }

	## make sure the output directory exists
    if (! -e "$options{'output_dir'}") {
    	$logger->logdie("The $options{'output_dir'} output directory passed could not be read or does not exist");
    }

}

# Subroutine to read files
sub read_file($) {
    my $filename = shift;
    my @lines;
    open(FH , "<$filename")  || $logger->logdie("Could not open $filename file for reading.$!");
    @lines = <FH>;
    close(FH);
    return(@lines);
}

sub _connect($$$$) {
    my ($user, $password, $db, $server) = @_;
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
