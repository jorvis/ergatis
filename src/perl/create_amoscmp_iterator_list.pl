#!/usr/bin/perl

eval 'exec /usr/bin/perl  -S $0 ${1+"$@"}'
    if 0; # not running under some shell
BEGIN{foreach (@INC) {s/\/usr\/local\/packages/\/local\/platform/}};
use lib (@INC,$ENV{"PERL_MOD_DIR"});
no lib "$ENV{PERL_MOD_DIR}/i686-linux";
no lib ".";

=head1 NAME

create_amoscmp_iterator_list.pl - Default output is a workflow iterator that can be used to iterator over input for amoscmp.

=head1 SYNOPSIS

USAGE: ./create_amoscmp_iterator_list.pl --sequence_list=/path/to/sequence/files/list --qual_list=/path/to/qual/files/list 
										 --reference_genome=/path/to/reference/genome --reference_genome_map=/path/to/reference/genome/map
										 --output=/path/to/output/iterator/xml

=head1 OPTIONS

B<--sequence_list, -s>
	One of two input files needed by toAmos; the list of sequence files that should be use to create an AMOS compatible afg format
	
B<--qual_list, -q>
	One of two input files needed by toAmos; the list of quality files that should be used to create an AMOS compatible afg file.	

B<--reference_genome, -r>
	A reference genome for reference based alignment.
	
B<--reference_map, -m>
	A reference genome map for aligning reference genomes against multiple inputs.		
	
B<--output, -o>
	Output iterator file.
	
B<--log, -l>
	Optional. Log file.

B<--debug, -d>
	Optional. Debug level.

B<--help, -h>
	Print perldocs for this script.
	
=head1 DESCRIPTION

Creates an ergatis/workflow iterator list file for a distributed AMOScmp job. This iterator contains the parameters needed to run toAmos,
converting sequence files (in FASTA format) and quality files to AMOScmp compatible afg format, and AmosCMP. 

=head1 INPUT

Two file lists, one containins sequence files (in FASTA format) and the other containing quality files

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
my $sequence_list = $options{'sequence_list'};
my $quality_list = $options{'qual_list'};
my $output_iterator = $options{'output'};
my $ref_genome = $options{'reference_genome'} if defined ($options{'reference_genome'});
my $ref_map_file = $options{'reference_map'} if defined ($options{'reference_map'});

## Parse sequences and quality file lists
my $sequences = &parse_file_list($sequence_list, ".fsa");
my $quals = &parse_file_list($quality_list, ".qual");
my $ref_map = &parse_reference_map($ref_map_file);

open (OUTFILE, "> $output_iterator") or $logger->logdie("Could not open output file $output_iterator: $!");
print OUTFILE '$;I_FILE_BASE$;' . "\t" .
			  '$;SEQUENCE_FILE$;' . "\t" .
			  '$;QUAL_FILE$;' . "\t" .
			  '$;REF_GENOME$;' . "\n";
			  
			  
foreach my $filename (keys %$sequences) {
	## Check if a quality file exists for this sequence file
	if ( defined($quals->{$filename}) ) {
		my $reference = get_reference_genome($ref_genome, $ref_map, $sequences->{$filename});
		## If we have a pair of files then write out the line to the iterator file.
		print OUTFILE "$filename\t$sequences->{$filename}\t$quals->{$filename}\t$reference\n";
	} else {
		$logger->logdie("Quality file does not exists for the sequence file $sequences->{$filename}");
	}
}

close (OUTFILE);

#########################################################################
#                                                                       #
#                           SUBROUTINES                                 #
#                                                                       #
#########################################################################

sub get_reference_genome {
	my ($ref, $ref_map, $filename) = @_;
	my $r_genome;
	
	if ( defined($ref) ) {
		$r_genome = $ref;	
	} else {
		$r_genome = $ref_map->{$filename};
	}
	
	unless ( defined ($r_genome) ) {
		$logger->logdie("Reference not found for file $filename");
	}
	
	return $r_genome;
}

sub parse_reference_map {
	my $ref_map_file = shift;
	my $ref_map = ();
	
	open (REFMAP, $ref_map_file) or $logger->logdie("Could not open reference map $ref_map_file: $!");
	while (my $line = <REFMAP>) {
		chomp ($line);
		my ($input, $reference) = split("\t", $line);
		
		## Check both of these files exist and are readble
		verify_file($input, $reference);
		
		unless ( exists($ref_map->{$input}) ) {
			$ref_map->{$input} = $reference;
		} else {
			$logger->logdie("Duplicate input files $input");
		}
	}
	
	return $ref_map;
}

## Parse the files out from the list passed in. Returns a hash
## key'd on filename.
sub parse_file_list {
	my ($list, $suffix) = @_;
	my $files = ();
	
	open (FLIST, $list) or $logger->logdie("Could not open file list $list: $!");
	while (my $line = <FLIST>) {
		chomp($line);
		
		if ( &verify_file($line) ) {
			my $filename = basename($line, $suffix);
			$files->{$filename} = $line;			
		}
	}
	
	close(FLIST);
	return $files;
}

## Verify that a file exists, is readable, and has content.
sub verify_file {
    my @files = @_;
    
    foreach my $file (@files) {
    	next if ( (-e $file) && (-r $file) && (-s $file) );
    	
    	if      (!-e $file) { $logger->logdie("File $file does not exist")   }   
        elsif   (!-r $file) { $logger->logdie("File $file is not readable")  }
    }
    
    return 1;
}

## Run system command and check for successful return value.
sub run_sys_cmd {
	my $cmd = shift;
	system($cmd);
	my $success = 1 if ($? >> 8 == 0);
	
	$logger->logdie("Could not execute command $cmd") unless ($success);
}

sub parse_options {
	my %opts = ();
	
	GetOptions (\%opts,
				 'sequence_list|s=s',
				 'qual_list|q=s',
				 'reference_genome|r=s',
				 'reference_map|m=s',
				 'output|o=s',
				 'log|l=s',
				 'debug|d=s',
				 'help|h') || pod2usage();
				 
	if ( $opts{'help'} ) {
	    pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDERR} );
	}
	
	## Set logger
	my $logfile = $opts{'log'} || Ergatis::Logger::get_default_filename();
	$logger = new Ergatis::Logger( 'LOG_FILE' 	=>	$opts{'log'},
        						   'LOG_LEVEL'	=>	$opts{'debug'} );
    $logger = Ergatis::Logger::get_logger();
    
    ## Check some of the parameters
    defined ($opts{'sequence_list'}) || $logger->logdie("Please specify a list of sequence files");
    defined ($opts{'qual_list'}) || $logger->logdie("Please specify a list of quality files");
    defined ($opts{'output'}) || $logger->logdie("Please specify a valid output directory");

	# We need something in input_file or input_list otherwise die and inform the user.
	$logger->logdie("Please provide either a single reference genome or a reference map") unless ( defined($opts{'reference_genome'}) || defined($opts{'reference_map'}) );
    
    return %opts;				 
}
				