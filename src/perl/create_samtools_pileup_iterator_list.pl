#!/usr/bin/perl

eval 'exec /usr/bin/perl  -S $0 ${1+"$@"}'
    if 0; # not running under some shell

eval 'exec /usr/bin/perl  -S $0 ${1+"$@"}'
    if 0; # not running under some shell
BEGIN{foreach (@INC) {s/\/usr\/local\/packages/\/local\/platform/}};
use lib (@INC,$ENV{"PERL_MOD_DIR"});
no lib "$ENV{PERL_MOD_DIR}/i686-linux";
no lib ".";

=head1 NAME

create_samtools_pileup_iterator_list.pl - Creates an iterator file list associated SAM files with a 
                                          corresponding reference FASTA file.
                                
=head1 SYNOPSIS

USAGE: ./create_samtools_pileup_iterator_list.pl --reference=/path/to/reference/file
                                                 --bam_file_list=/path/to/bam/file/list
                                                 --output=/path/to/output/iterator/file
=head1 OPTIONS

B<--reference, -r>
    A file or file list containing a single reference FASTA file.
    
B<--bam_file_list, -s> 
    A file or file list containing an input BAM file.

B<--output, -o>
    Desired output iterator file list.

B<--help>
    Print perldocs for this script.
    
=head1 DESCRIPTION

Creates an ergatis/workflow iterator list file for a distributed sam_pileup job. The iterator list file
requires an input SAM file coupled to a single reference FASTA file.

=head1 INPUT

A single SAM file or list of SAM files alongside with a single FASTA reference sequence.

=head1 OUTPUT
 
An ergatis iterator list

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
my $references = $options{'reference'};
my $bam_file_list = $options{'bam_file_list'};
my $output_iterator = $options{'output'};

my $reference_file = parse_reference_file($references);
my $bam_files = parse_bam_file_list($bam_file_list);
write_output_iterator($reference_file, $bam_files, $output_iterator);

#########################################################################
#                                                                       #
#                           SUBROUTINES                                 #
#                                                                       #
#########################################################################
 
#--------------------------------------------------------------------- 
# Creates an Ergatis output iterator in the following format:
#  
# $;I_FILEBASE$;\t$;I_FILE_NAME$;\t$;I_FILE_PATH$;\t$;REF_SEQ$;
# sam_file_1\tsam_file_1.sam\t/tmp/sam_file_1.sam\tref_fasta.fsa
#
#---------------------------------------------------------------------
sub write_output_iterator {
    my ($ref, $bams, $iter) = @_;

    open (OUTITER, "> $iter") or $logger->logdie("Could not open output iterator $iter for writing: $!");
    print OUTITER '$;I_FILE_BASE$;' . "\t" .
                  '$;I_FILE_NAME$;' . "\t" .
                  '$;I_FILE_PATH$;' . "\t" .
                  '$;REF_FILE$;' . "\n";

    foreach my $bam_file (@{ $bams }) {
        my $filename = basename($bam_file);
        my $file_base = fileparse($bam_file, '\.(.*)');

        print OUTITER "$file_base\t$filename\t$bam_file\t$ref\n";
    }

    close OUTITER;
}

#--------------------------------------------------------------------- 
# Parses the BAM file list or returns the single BAM file that was
# passed into this component.
#---------------------------------------------------------------------
sub parse_bam_file_list {
    my $bams = shift;
    my $bam_files;

    foreach my $bam_file (@{ $bams }) {
        open (BAMFILE, $bam_file) or $logger->logdie("Could not open BAM file $bam_file: $!");
        my @files = <BAMFILE>;
        chomp(@files);

        if (-e $files[0]) {
            push(@{ $bam_files }, @files); 
        } else {
            push(@{ $bam_files }, $bam_file);
        }
    }

    close (BAMFILE);
    return $bam_files;
}

#--------------------------------------------------------------------- 
# Parses the reference file or reference file list passed into the 
# component. If a file list is passed in verification is done to 
# ensure that only one file is present in the list.
#---------------------------------------------------------------------
sub parse_reference_file {
    my $refs = shift;
    my $ref_seq = undef;

    open (REFFILE, $refs->[0]) or $logger->logdie("Could not open reference $refs: $!");
    my @lines = <REFFILE>;
    
    if ($lines[0] =~ /^>/) {
        $ref_seq = $refs;
    } else {
        $logger->logdie("Only one reference FASTA can be passed into this component.") if (scalar @lines > 1);
        chomp($lines[0]);
        my $file = $lines[0];

        if (-e $file) {
            $ref_seq = $file;
        } else {
            $logger->logdie("File $file does not exist");
        }
    }

    close (REFFILE);
    return $ref_seq;
}

#--------------------------------------------------------------------- 
# Parse command-line arguments                                       
#---------------------------------------------------------------------
sub parse_options {
    my %opts = ();

    GetOptions(\%opts,
                'reference|r=s@',
                'bam_file_list|s=s@',
                'output|o=s',
                'help' ) || pod2usage();

    if ($opts{'help'}) {
        pod2usage ( { -exitval => 0, -verbose => 2, -output => \*STDERR } );
    }
    
    my $logfile = Ergatis::Logger::get_default_logfilename();
    my $debug = 4;
    $logger = new Ergatis::Logger( 'LOG_FILE'   =>  $logfile,
                                   'LOG_LEVEL'  =>  $debug );
    $logger = Ergatis::Logger::get_logger();

    defined ($opts{'reference'}) || $logger->logdie("Please specify a reference file.");
    defined ($opts{'bam_file_list'}) || $logger->logdie("Please specify a single BAM file or list of BAM files.");
    defined ($opts{'output'}) || $logger->logdie("Please specify the desired output iterator file.");
    
    # We don't know if we will receive our reference as a file or file list we
    # need to do a quick check here to ensure our array is not larger than two
    if (scalar @{ $opts{'reference'} } > 1) {
        $logger->logdie("Only one reference FASTA may be passed into this component");
    }

    return %opts
}
                                        
