#!/usr/bin/perl

=head1 NAME

process_crossbow_reference_genome.pl - Processes one or more reference genomes to be used in short read mapping via crossbow.

=head1 USAGE

./process_crossbow_reference_genome.pl 
        --bowtie_build_exec=/path/to/bowtie/directory
        --reference_genomes=/path/to/reference/genome(s)
        --output_dir=/path/to/desired/output/directory
       [--log=/desired/path/to/log/file
        --debug=/desired/log/level
        --help]
        
=head1 SYNOPSIS

B<--bowtie_build_exec, -b>
    The path to the bowtie-build executable.
    
B<--reference_genomes, -r>
    Either a single reference genome or a list of reference genomes separate by comma.
    
B<--output_dir, -o>
    The directory to place index files produced from bowtie-build
    
B<--log, -l>
    Optional. A log file to house any debug/warning/error statements.
    
B<--debug, -d>
    Optional. Level of logging to use. Can be set to between 1 - 4 (increasing in verbosity)
    
B<--help, -h>
    Print help documentation
    
=head1 DESCRIPTION

Script processes any reference genomes that will be used in short read mapping with Crossbow.
Reference genomes must be processed and indexed using the bowtie-build executable.

=head1 INPUT

Input is one or more reference genomes in FASTA file. These genomes must be processed in the
following order:

        1.) Each subsequent genome must be renamed to the following format chrX 
            where X = 1. 
        2.) Each reference genome should be a single sequence and the header line
            should be set to the following ">X" where X matches the X in the filename.
            
=head1 OUTPUT

A set of index files produced by bowtie.

=head1 CONTACT

    Cesar Arze
    carze@som.umaryland.edu
    
=cut

use strict;
use warnings;
use Pod::Usage;
use FileHandle;
use Log::Log4perl qw(:easy);
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);

#---------------------------------------------
# GLOBALS/COMMAND-LINE ARGUMENTS
#---------------------------------------------
my $logger;

my %options = &parse_options();
my $BOWTIE_BUILD_EXEC = $options{'bowtie_build_exec'};
my $references = $options{'reference_genomes'};
my $output_directory = $options{'output_dir'};

my @reference_genomes = &process_reference_genomes($references, $output_directory);

## Processed reference genomes must be run though the bowtie-build executable.
my $genome_string = join(",", @reference_genomes);
my $cmd = "$BOWTIE_BUILD_EXEC $genome_string $output_directory/index";
run_system_cmd($cmd);

###############################################################################
#####                          SUBROUTINES                                #####
###############################################################################

#---------------------------------------------
# process reference genome(s)
#---------------------------------------------
sub process_reference_genomes {
    my ($ref_genomes, $output_dir) = @_;
    my @processed_genomes;
    
    ## One or more genomes should be provided in a comma-delimited string
    my @raw_genomes = split(/,/, $ref_genomes);
    my $ref_index = 0;
    
    foreach my $raw_genome (@raw_genomes) {
        ## Check if this is a valid file first
        validate_file($raw_genome);
        
        ## Verify that the FASTA file only contains one sequence
        my $in_fh = FileHandle->new($raw_genome, "r") or $logger->logdie("Could not open reference genome $raw_genome: $!");
        my $seq_count = &sequence_count($in_fh);
        if ( $seq_count > 1 || $seq_count == 0) { $logger->logdie("Reference genome file $raw_genome contains more than one sequence.") }
        
        ## The header of the reference genome must be set to the numeric in the "chr" filename
        my $processed_genome = $output_dir . "/" . "chr" . $ref_index . ".fa";
	$logger->logdie("File $processed_genome already exists.") if (-e $processed_genome);
        open (GENOMEOUT, ">" . $processed_genome) or $logger->logdie("Could not write to processed reference genome file $processed_genome: $!");
        while (my $line = $in_fh->getline) {
            if ($line =~ /^>/) {
                print GENOMEOUT ">$ref_index\n";
            } else {
                print GENOMEOUT $line;
            }
        }
        
        close (GENOMEOUT);
        push (@processed_genomes, $processed_genome);
        $ref_index++;            
    }
    
    return @processed_genomes;
}

#---------------------------------------------
# count number of sequences in a FASTA file
#---------------------------------------------
sub sequence_count {
    my $in_ref = shift;
    my $seq_count = 0;
    while (my $line = $in_ref->getline) { $seq_count++ if ($line =~ /^>/) }
    seek ($in_ref,0, 0);
    return $seq_count;
}

#---------------------------------------------
# verify that a file exists and is readble
#---------------------------------------------
sub validate_file {
    my $file = shift;
 
    if      (!-e $file) { $logger->logdie("File $file does not exist")   }   
    elsif   (!-r $file) { $logger->logdie("File $file is not readable")  }
    elsif   (!-s $file) { $logger->logdie("File $file has zero content") }  
}

#-----------------------------------
# run a unix system command
#-----------------------------------
sub run_system_cmd {
    my $cmd = shift;
    my $res = `$cmd`;
    chomp($res);
    my $success = $? >> 8;

    unless ($success == 0) {
        $logger->logdie("Command \"$cmd\" failed: $res");
    }
}

#---------------------------------------------
# parse command-line arguments
#---------------------------------------------
sub parse_options {
    my %opts = ();
    GetOptions(\%opts, 
                'bowtie_build_exec|b=s',
                'reference_genomes|r=s',
                'output_dir|o=s',
                'log|l:s',
                'debug|d:s',
                'help') || pod2usage();
    
	&pod2usage( {-exitval => 1, -verbose => 1, -output => \*STDOUT}) if ($opts{'help'} );			

	## Configure Log4perl
	my $debug_lvl = $opts{'debug'} ||= "DEBUG";
	my $log_file = $opts{'log'} ||= "/tmp/demultiplex_sff_with_read_lists.pl.log";
	Log::Log4perl->easy_init( { level => $debug_lvl, file => $log_file } );
	$logger = get_logger();        
	
	## Check that necessary parameters are defined
	defined ($opts{'bowtie_build_exec'}) || $logger->logdie("Please provide the path to a bowtie-build executable");
	defined ($opts{'reference_genomes'}) || $logger->logdie("Please provide or more reference genomes in FASTA format");
	defined ($opts{'output_dir'}) || $logger->logdie("Please provide an output directory.");
    
    return %opts;
}
