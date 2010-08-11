#!/usr/bin/perl

=head1 NAME

transform_WWARN_input.pl - Transforms input data to WWARN to a common format.

=head1 SYNOPSIS

./tag_data.pl 
        --input_file=/path/to/input/file
        --input_file_type=<input file type>
        --tag-name=<clovr tag name>
       [--log=/path/to/log/file
        --debug=/debug/lvl
        --help]
        
=head1 PARAMETERS

B<--input_file, -i>
	Input files which should be tagged by CloVR. Input files can be either individual 
	files, a directory containing many files or a file containing a list of files
    
B<--input_file_type, -t>
	The type of file specified in the input_file parameter. Valid options here are DIR, FILE,
	and LIST
	
B<--tag-name, -n>
	The unique name to tag these data with in CloVR.	

B<--log, -l>
    Optional. Log file.
    
B<--debug, -d>
    Optional. Debug level.
    
B<--help>
    Prints this documentaiton.
    
=head1 DESCRIPTION            

A wrapper script that facilitates the execution of the CloVR tagData.py script through ergatis.
This wrapper script is meant to allow for custom tagging of output data from a pipeline.

=head1 INPUT

Input is a single input file, either a directory, file, or file list that contains data to be tagged.

=head1 CONTACT

    Cesar Arze
    carze@som.umaryland.edu

=cut

use strict;
use warnings;
use strict;
use warnings;
use Pod::Usage;		
use Ergatis::Logger;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);

#----------------------------------------------------------
# GLOBALS/COMMAND-LINE OPTIONS
#----------------------------------------------------------
my $logger;

my %options = parse_options();
my $input = $options{'input_file'};
my $input_file_type = $options{'input_file_type'};
my $tag_name = $options{'tag-name'};

my $cmd = qq{/opt/vappio-py/vappio/cli/tagData.py --name=local --tag-name="$tag_name" --overwrite --recursive };

if ( uc($input_file_type) eq "LIST" ) {
	open (INFILE, $input);
	chomp( my @files = <INFILE> );
	close (INFILE);

	$cmd .= join(" ", @files);
} else {
	$cmd .= "$input";
}

run_system_cmd($cmd);

###############################################################################
#####                          SUBROUTINES                                #####
###############################################################################

sub run_system_cmd {
    my $cmd = shift;
   
    system($cmd);
    my $success = $? >> 8;

    unless ($success == 0) {
        $logger->logdie("Command \"$cmd\" failed.");
    }   
}

sub parse_options {
	my %opts = ();
	GetOptions(\%opts,
				'input_file|i=s',
				'input_file_type|t=s',
				'tag-name|n=s',
                'log|l:s',
                'debug|d:s',
                'help') || pod2usage();
                 
	&pod2usage( {-exitval => 1, -verbose => 1, -output => \*STDOUT}) if ($opts{'help'} );		
	
    ## Configure logger
    my $logfile = $opts{'log'} || Ergatis::Logger::get_default_logfilename();
    $logger = new Ergatis::Logger( 'LOG_FILE'   =>  $opts{'log'},
                                   'LOG_LEVEL'  =>  $opts{'debug'} );
    $logger = Ergatis::Logger::get_logger();
    
    ## Make sure our parameter are declared correctly
    defined($opts{'input_file'}) || $logger->logdie('Please specify a valid input file');
    defined($opts{'input_file_type'}) || $logger->logdie('Please specify a valid input file type [DIR, FILE, FILE LIST]');
    defined($opts{'tag-name'}) || $logger->logdie('Please specify a valid tag name');
	
	return %opts;
}
