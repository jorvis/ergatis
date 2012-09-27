#!/usr/bin/perl -w

=head1 NAME

 predict_operons.pl - wrapper script for predict_operons.sh

=head1 SYNOPSIS

 perl predict_operons.pl -b /bin/dir -o /output/dir -w /wig/dir/ -g /path/to/file.gff3 -r /path/to/reference.fa -s sample_name

=head1 OPTIONS

	--bin_dir|b	-	The path of the bin directory where the predict_transcripts.sh, find_transcripts.pl, and find_intergenic_background_cutoff.py are stored
	--output_dir|o	-	The output directory to store files
	--wig_dir|w	-	The directory where WIG files are stored (typically will be the output directory of the bam2bigwig.visualization Ergatis component
	--gff|g		-	The path to the GFF3 file
	--reference|r	-	A reference fasta file path
	--sample|s	-	The name of the sample (typically will be the abbreviation used in Ergatis ProkPipe)

=head1 DESCRIPTION

 This is a wrapper script that will provide these arguments as parameters to predict_operons.sh.  That shell script in turn, calls several other scripts.
 Fastalen.awk will create a tab-separated (.tsv) file of contigs and contig lengths from the reference fasta file
 Find_intergenic_background_cutoff.py will produce a depth-of-coverage cut-off intended for transcript finding
 Find_transcripts.pl will take the cutoff output file and find transcripts from an alignment based on coverage islands

=head1 OUTPUT

 find_transcripts.pl will output a GFF3 file which will list the coordinates of all transcripts, and another GFF3 which will list only those that are not known (do not overlap with known)

=head1 CONTACT

	Shaun Adkins
	sadkins@som.umaryland.edu

=cut

use strict;
use warnings;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Pod::Usage;


my %options;
GetOptions(\%options,
	   'bin_dir|b=s',
	   'output_dir|o=s',
	   'wig_dir|w=s',
	   'gff|g=s',
	   'reference|r=s',
	   'sample|s=s', 
	   'help|h'
	  ) or pod2usage();

pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDERR} ) if ($options{'help'});

&check_parameters(%options);

my $shell_script = $options{'bin_dir'} . '/predict_transcripts.sh';

`sh $shell_script $options{'bin_dir'} $options{'output_dir'} $options{'wig_dir'} $options{'reference'} $options{'gff'} $options{'sample'}`;

sub check_parameters {
    my $options = shift;
    # make sure required arguments were passed
    my @required = qw( bin_dir output_dir wig_dir gff reference sample );
    for my $option ( @required ) {
        unless  ( defined $options{$option} ) {
            die "--$option is a required option";
        }
    }
}
