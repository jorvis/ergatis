#!/usr/bin/perl
=head1  NAME 

create_update_euk_model_sequences_iterator_list.pl - Default output is a workflow iterator that
can be used to iterator over a set of database-asmbl_id values.  Can read in the standard legacy2bsml
control file.

=head1 SYNOPSIS

USAGE:  create_update_euk_model_sequences_iterator_list.pl --control_file=/tmp/afu1.control-file.dat --output=/tmp/outfile.txt

=head1 OPTIONS

=item *

B<--debug,-d> Debug level.  Use a large number to turn on verbose debugging. 

=item *

B<--log,-l> Log file

=item *

B<--help,-h> This help message

=head1   DESCRIPTION - TBA

Need to provide some description here

=head1   CONTACT

Jay Sundaram

sundaram@jcvi.org

=head1   FYI

I've recycled portions of the create_legacy2bsml_iterator_list.pl script.

=cut

no lib "$ENV{PERL_MOD_DIR}/i686-linux";
no lib ".";

use strict;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use File::Basename;

BEGIN {
    require '/usr/local/devel/ANNOTATION/cas/lib/site_perl/5.8.5/Workflow/Logger.pm';
    import Workflow::Logger;
}


umask(0000);

my %options = ();

my $results = GetOptions (\%options, 
                          'control_file|f=s', 
			  'output|o=s',
			  'log|l=s',
                          'debug=s', 
                          'help|h' ) || pod2usage();

my $logfile = $options{'log'} || Workflow::Logger::get_default_logfilename();
my $logger = new Workflow::Logger('LOG_FILE'=>$logfile,
				  'LOG_LEVEL'=>$options{'debug'});
$logger = Workflow::Logger::get_logger();

## Check critical command-line arguments

if (! &verifyControlFile($options{'control_file'})){
    $logger->logdie("User did not specify the --control_file argument");
}

if (! defined($options{'output'})){
    $logger->logdie("User did not specify the --output argument");
}

## Declare master array
my $iteratorconf = [];

## Read in the information from asmbl_file OR asmbl_list
&getListFromFile($iteratorconf,$options{'control_file'});

## Output the lists
&createOutputListFile($iteratorconf, $options{'output'});

## Exit trivialities
print "$0 script execution completed\n";
print "The log file is '$logfile'\n";
exit(0);
						     
#---------------------------------------------------------------------------------------------------------
#
#                           END OF MAIN  --  SUBROUTINES FOLLOW
#
#---------------------------------------------------------------------------------------------------------

=over 4

=item getListFromFile()

B<Description:> 

B<Parameters:> $iteratorconf (array reference), $file (scalar)

B<Returns:> None

=back

=cut

sub getListFromFile {

    my ($iteratorconf, $file) = @_;

    my $contents = &getFileContents($file);

    my $dataLookup = &getDataLookup($contents, $file);

    foreach my $database (sort keys %{$dataLookup} ){ 

	foreach my $asmbl_id (sort {$a <=> $b } @{$dataLookup->{$database}}){
	    
	    my $uniqueId = $database . '_' . $asmbl_id;
	    
	    
	    my @tempArray = ( $uniqueId, 
			      $database, 
			      $asmbl_id
			      );
	    
	    push(@{$iteratorconf}, \@tempArray);
	    
	}
    }
}


=over 4

=item createOutputListFile()

B<Description:> 

B<Parameters:> $iteratorconf (array reference), $file (scalar)

B<Returns:> None

=back

=cut

sub createOutputListFile {

    my ($iteratorconf, $output_file) = @_;

    open (OUTPUT_FILE, "+>$output_file") || $logger->logdie("Could not open output file '$output_file' in output mode:$!");
    
    print OUTPUT_FILE '$;UNIQUE_KEY$;' . "\t".
    '$;DATABASE$;' . "\t".
    '$;ASMBL_ID$;' . "\n";

    foreach my $arrayRef (@{$iteratorconf}){
	
	print OUTPUT_FILE join("\t",@{$arrayRef}),"\n";
	
    }

    close FILE;

}

=over 4

=item verifyControlFile()

B<Description:> 

B<Parameters:> $file (scalar)

B<Returns:> 1 - true (scalar) OR 0 - false (scalar)

=back

=cut

sub verifyControlFile {

    my $file = shift;

    if ((defined($file)) && (-e $file) && (-r $file) && (-s $file)) {
	# control file was defined, exists, has read permissions
	# and has content
	return 1;
    }
    else {
	if (!defined($file)){
	    $logger->logdie("--control_file was not defined");
	}
	if (!-e $file){
	    $logger->logdie("control file '$file' does not exist");
	}
	if (!-r $file){
	    $logger->logdie("control file '$file' does not have read permissions");
	}
	if (!-s $file){
	    $logger->logdie("control file '$file' has zero content");
	}
    }
}

=over 4

=item getFileContents()

B<Description:> 

B<Parameters:> $file (scalar)

B<Returns:> reference to array

=back

=cut

sub getFileContents {

    my $file = shift;

    open (CONTROLFILE,  "<$file") or $logger->logdie("Could not open file '$file' in read mode: $!");
    
    my @lines = <CONTROLFILE>;

    chomp @lines;

    return \@lines;
}


=over 4

=item getDataLookup()

B<Description:> 

B<Parameters:> $file (scalar)

B<Returns:> reference to array

=back

=cut

sub getDataLookup {

    my ($contents, $control_file) = @_;

    my $lookup = {};

    ## Should only process unique sets of database-asmbl_id values
    my $uniqLookup = {};

    ## Keep track of number of unique values encountered
    my $uniqValCtr=0;

    ## Keep track of duplicate values
    my $dupLookup = {};

    ## Keep track of number of duplicate values encountered
    my $dupCtr=0;


    ## Keep track of the number of control file lines get processed
    my $linectr=0;

    my $database;

    foreach my $line (@{$contents}){

	$linectr++;

	if ($line =~ /^\s*$/){
	    next; # skip blank lines
	}
	elsif ($line =~ /^\#/){
	    next; # skip comment lines
	}
	elsif ($line =~ /^\-\-/){
	    next; # skip -- lines
	}
	else{

	    if ($line =~ /^database:(\S+)/){
		## This is the line which identifies the annotation database to be processed
		$database = $1;
	    }
	    elsif ($line =~ /^\s*(\d+)\s*$/){
		## All other lines should contain asmbl_id values
		my $asmbl_id = $1;
		
		my $key = $database . '_' . $asmbl_id;
		
		if (! exists $uniqLookup->{$key}){
		    push(@{$lookup->{$database}}, $asmbl_id);
		    $uniqValCtr++;
		}
		else {
		    $dupLookup->{$key}++;
		    $dupCtr++;
		}
		
		$uniqLookup->{$key}++;
	    }
	    else {
		$logger->logdie("Could not parse line number '$linectr' - line was '$line'");
	    }
	}
    }

    if ($dupCtr>0){
	$logger->warn("Encountered a number of repeated values in control file '$control_file':");
	foreach my $key (sort keys %{$dupLookup}){
	    my ($database, $asmbl_id) = split("_", $key);
	    $logger->warn("database '$database' asmbl_id '$asmbl_id' was encountered '$dupLookup->{$key}' times");
	}
    }

    if ($uniqValCtr == 0){
	$logger->logdie("Did not encounter any unique values in the control file '$control_file':");
    }

    return $lookup;
}
