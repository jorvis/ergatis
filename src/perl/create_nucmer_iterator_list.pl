#!/usr/bin/perl

eval 'exec /usr/bin/perl  -S $0 ${1+"$@"}'
    if 0; # not running under some shell
BEGIN{foreach (@INC) {s/\/usr\/local\/packages/\/local\/platform/}};
use lib (@INC,$ENV{"PERL_MOD_DIR"});
no lib "$ENV{PERL_MOD_DIR}/i686-linux";
no lib ".";

=head1 NAME

create_nucmer_iterator_list.pl - Default output is a workflow iterator that can be used to iterate over input for nucmer

=head1 SYNOPSIS

USAGE: ./create_nucmer_iterator_list.pl --sequence_list=/path/to/sequence/file/list --reference_genome=/path/to/reference/genome
										--reference_genome_list=/path/to/reference/genome/list --output=/path/to/output/iterator
										
=head1 OPTIONS

B<--sequence_file, -i>
    A single sequence file.

B<--sequence_list, -s>
	A list of sequence files that should be aligned using nucmer.
	
B<--reference_genome, -r>
	A reference genome to align the target sequences against.
	
B<--reference_genome_list, -rl>
	A list of reference genome that every target input sequence should be aligned against.

B<--output, -o>
	Output iterator file.
	
B<--log, -l>
	Optional. Log file.
	
B<--debug, -d>	
	Optional. Debug level.
	
B<--help, -h>
	Print perldocs for this script.
	
=head1 DESCRIPTION

Creates an ergatis/workflow iterator list file for a distributed nucmer job. This iterator contains the parameters needed to run nucmer.

=head1 INPUT

A list of sequence files in FASTA format and a single reference genome (in FASTA format) AND/OR a list of reference genomes.

=head1 OUTPUT

The output ergatis iterator list will be written to the file specified by the --output parameter.

=head1 CONTACT

	Cesar Arze
	carze@som.umaryland.edu

=cut

use strict;
use warnings;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Pod::Usage;
use File::Basename;
use Ergatis::Logger;
umask(0000);
my $logger;

my %options = &parse_options();
my $sequence_file = $options{'sequence_file'};
my $sequence_list = $options{'sequence_list'};
my $ref_genome = $options{'reference_genome'} if defined ($options{'reference_genome'});
my $ref_genome_list = $options{'reference_genome_list'} if defined ($options{'reference_genome_list'});

my @sequences = &parse_sequence_list($sequence_list, $sequence_file);
my @references = &parse_reference_genomes($ref_genome, $ref_genome_list);

open (OUTFILE, "> $options{'output'}") or $logger->logdie("Could not open output iterator $options{'ouput'} for writing: $!");
print OUTFILE '$;UNIQUE_BASE$;' . "\t" .
			  '$;REFERENCE_GENOME$;' . "\t" .
              '$;I_FILE_BASE$;' . "\t" .
			  '$;I_FILE_PATH$;' . "\n"; 

foreach my $sequence (@sequences) {
	my $prefix = fileparse($sequence, '.fsa|.fna|.fasta');
	foreach my $ref_genome (@references) {
        my $ref_filename = fileparse($ref_genome, '\.[^\.]*');
        $ref_filename =~ s/\.nuc//; # Hack for oral metagenomics pipeline
		print OUTFILE "$prefix\.$ref_filename\t$ref_genome\t$prefix\t$sequence\n";
	}
}						  

close (OUTFILE);

#########################################################################
#                                                                       #
#                           SUBROUTINES                                 #
#                                                                       #
#########################################################################

## Parse the files out from the list passed in. Returns a hash
## key'd on filename.
sub parse_sequence_list {
	my ($list, $file) = @_;
	my @files = ();
	
    ## Handle a single FASTA file being passed in...
    push (@files, $file) if ( defined($file) && &verify_file($file) );

    if ( defined($list) ) {
    	open (FLIST, $list) or $logger->logdie("Could not open file list $list: $!");
    	while (my $line = <FLIST>) {
    		chomp($line);
		
    		if ( &verify_file($line) ) {
    			push (@files, $line);			
    		}
    	}
    }
	
	close(FLIST);
	return @files;
}

## Verify that a file exists, is readable, and has content.
sub verify_file {
    my @files = @_;
    
    foreach my $file (@files) {
    	next if ( (-e $file) && (-r $file) && (-s $file) );
    	
    	if      (!-e $file) { $logger->logdie("File $file does not exist")   }   
        elsif   (!-r $file) { $logger->logdie("File $file is not readable")  }
        elsif   (!-s $file) { $logger->logdie("File $file has zero content") }
    }
    
    return 1;
}

## Handle either a single reference genome our a list of reference genomes.
## If a list is passed in we will iterate our input over every reference in
## the list.
sub parse_reference_genomes {
	my ($ref_genome, $ref_list) = @_;
	my @refs = ();
	
	push (@refs, $ref_genome) if ( defined ($ref_genome) && &verify_file($ref_genome) );
	if ( &verify_file($ref_list) ) {
		open(REFLIST, $ref_list) or $logger->logdie("Could not open reference genome list $ref_list: $!");
		
		while (my $line = <REFLIST>) {
			chomp($line);
			push(@refs, $line) if ( &verify_file($line) );
		}
		
		close (REFLIST);
	}
	
    $logger->logdie("No reference genomes provided.") if (scalar @refs == 0);
	return @refs;
}

sub parse_options {
	my %opts = ();
	
	GetOptions (\%opts,
                'sequence_file|i=s',
				'sequence_list|s=s',
				'reference_genome|r=s',
				'reference_genome_list|rl=s',
				'output|o=s',
				'log|l=s',
				'debug|d=s',
				'help') || pod2usage();
				
	if ( $opts{'help'} ) {
	    pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDERR} );
	}
	
	## Set logger
	my $logfile = $opts{'log'} || Ergatis::Logger::get_default_filename();
	$logger = new Ergatis::Logger( 'LOG_FILE' 	=>	$opts{'log'},
        						   'LOG_LEVEL'	=>	$opts{'debug'} );
    $logger = Ergatis::Logger::get_logger();
    
    ## Check certain parameters to make sure they are defined...
    # defined ($opts{'sequence_list'}) || $logger->logdie("Please specify a list of sequence files");
    defined ($opts{'output'}) || $logger->logdie("Please specify a valid output file");

	# We need something in input_file or input_list otherwise die and inform the user.
    $logger->logdie("Please provide either a single input filr or a list of input files") unless ( defined($opts{'sequence_file'}) || defined($opts{'sequence_list'}) );

    # We need something in reference_genome or reference_genome_list otherwise die and inform the user.
	$logger->logdie("Please provide either a single reference genome or a list of reference genomes") unless ( defined($opts{'reference_genome'}) || defined($opts{'reference_genome_list'}) );
			
	return %opts;			
}											
