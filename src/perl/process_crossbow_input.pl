#!/usr/bin/perl

=head1 NAME

process_crossbow_input.pl - Processes crossbow input

=head1 USAGE

./process_crossbow_input.pl 
        --input_file=/path/to/input/text/file
        --output_file=/desired/path/to/output/manifest/file
       [--log=/desired/path/to/log/file
        --debug=/debug/level
        --help]
        
=head1 SYNOPSIS

B<--input_file, -i>
    A tab-delimited text file containing the path to one or more FASTA files.
    
B<--output_file, -o>            
    Output manifest file for use in crossbow pipeline.
    
B<--log, -l>
    Optional. Log file to capture debugging statements.
    
B<--debug, -d>
    Optional. Debug level.
    
B<--help>
    Print out help documentation
    
=head1 DESCRIPTION

Creates a Crossbow manifest file.

=head1 INPUT

A tab-delimited list text file containing the path to one or more input FASTA files. If paired reads exist 
they should be separated by tabs:

/path/to/paired_read1.fasta\t/path/to/paired_read2.fasta

=head1 OUTPUT

Crossbow manifest file.

=head1 CONTACT

    Cesar Arze
    carze@som.umaryland.edu
                        
=cut

use strict;
use warnings;
use Pod::Usage;
use File::Basename;
use FileHandle;
use Log::Log4perl qw(:easy);
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);

#---------------------------------------------
# GLOBALS/COMMAND-LINE ARGUMENTS
#---------------------------------------------
my $logger;

my %options = &parse_options();
my $input_file = $options{'input_file'};
my $output_manifest = $options{'output_file'};
my $project_name = $options{'project_name'};

## Pull down and process hostname
my $hostname = `hostname`;
chomp($hostname);
$hostname = "hdfs://" . $hostname;

open (INPUT, $input_file) or $logger->logdie("Cannot open input file $input_file: $!");
open (OUTPUT, ">" . $output_manifest) or $logger->logdie("Cannot write to manifest file $output_manifest: $!");

while (my $line = <INPUT>) {
    chomp($line);
    my @reads = split(/\t/, $line);
    
    foreach my $read (@reads) {
        ## Make sure file is exists and is well-formed
        validate_file($read); 
        my $filename = basename($read);
        
        ## Upload the read to hadoop
        system("hadoop fs -put $read /users/clovr/$project_name/reads");
        
        ## Print out to manifest file
        print OUTPUT $hostname . "/users/clovr/" . $project_name . "/reads/" . $filename . " 0 ";
    }
    
    print OUTPUT "\n";
}

close (INPUT);
close (OUTPUT);

###############################################################################
#####                          SUBROUTINES                                #####
###############################################################################

#---------------------------------------------
# verify that a file exists and is readble
#---------------------------------------------
sub validate_file {
    my $file = shift;
 
    if      (!-e $file) { $logger->logdie("File $file does not exist")   }   
    elsif   (!-r $file) { $logger->logdie("File $file is not readable")  }
    elsif   (!-s $file) { $logger->logdie("File $file has zero content") }  
}


sub parse_options {
    my %opts = ();
    GetOptions(\%opts,
                'input_file|i=s',
                'output_file|o=s',
                'project_name|p=s',
                'log|l:s',
                'debug|d:s',
                'help') || pod2usage();
                
    &pod2usage( {-exitval => 1, -verbose => 1, -output => \*STDOUT}) if ($opts{'help'} );			

	## Configure Log4perl
	my $debug_lvl = $opts{'debug'} ||= "DEBUG";
	my $log_file = $opts{'log'} ||= "/tmp/create_crossbow_manifest.pl.log";
	Log::Log4perl->easy_init( { level => $debug_lvl, file => $log_file } );
	$logger = get_logger();        
	
	## Check to verify mandatory parameters are defined
	defined ($opts{'input_file'}) || $logger->logdie("Please provide a valid input file.");
	defined ($opts{'output_file'}) || $logger->logdie("Please provide a path to a desired output manifest file.");
	defined ($opts{'project_name'}) || $logger->logdie("Please provide a valid project name.");
	
	return %opts;
}                        
