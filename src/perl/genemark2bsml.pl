#!/usr/bin/perl

=head1 EXPECTED INPUT

    Predicted genes/exons

    Gene Exon Strand Exon           Exon Range     Exon      Start/End
      #    #         Type                         Length       Frame

      1     1   +  Terminal      6168      6449     282          1 3

      2     3   -  Terminal     13450     13528      79          3 3
      2     2   -  Internal     16097     16311     215          1 2
      2     1   -  Initial      16436     16468      33          1 3

      3     1   +  Initial      19541     19632      92          1 2
      3     2   +  Terminal     19755     20169     415          3 3

      4     1   +  Initial      34531     34622      92          1 2
      4     2   +  Internal     34745     34967     223          3 3
      4     3   +  Terminal     35854     35982     129          1 3

=cut

use strict;
use warnings;

use Getopt::Long qw(:config no_ignore_case bundling);
use File::Basename;
use IO::File;
use GenePredictionUtils::GeneMarkBsmlGenerator;

use constant
{
	EXON_LABELS	=>	{	"Initial"	=>	1,
					"Internal"	=>	1,
					"Terminal"	=>	1,
					"Single"	=>	1
				}
};

my $in			= new IO::File->fdopen(fileno(STDIN), "r");
my $out			= "/dev/stdout";
my $fasta		= undef;
my $project		= "unknown";
my $gene_finder_name	= "GeneMark.hmm";
my $sequence_id = undef;

&parse_options;
&create_bsml;

sub print_usage
{
	my $progname = basename($0);
	die << "END";
usage: $progname [-input|i <input_genscan_data>] [-output|o <output_bsml>]
	[-fasta_file|f <fasta_data>] [-project|p <project_name>] [-sequence_id|s <sequence_id>]
	[-help|h]
END
}

sub parse_options
{
	my %opts = ();
	GetOptions(\%opts, "input|i=s", "output|o=s", "fasta_file|f=s",
		   "project|p=s", "help|h", "sequence_id|s=s");
	&print_usage if $opts{help};
	$in->open($opts{input}, "r") or
		die "Error reading genscan input $opts{input}: $!"
		if $opts{input};
	$out = $opts{output} if $opts{output};
	$fasta = $opts{fasta_file} if $opts{fasta_file};
	$project = $opts{project} if $opts{project};
    $sequence_id = $opts{sequence_id} if $opts{sequence_id};
}

sub create_bsml
{
	my $bsml_generator = new GenePredictionUtils::GeneMarkBsmlGenerator;
	my %gene_data = ();
	while (my $line = <$in>) {
		chomp $line;
		next if $line =~ /^Sequence name/;
		next if !length($line);
		$line =~ s/^\s+//g;
		my @tokens = split /\s+/, $line;
		next if $tokens[0] !~ /^(\d+)/;
		if (exists EXON_LABELS->{$tokens[3]}) {
			my $from = $tokens[4] - 1;
			my $to = $tokens[5];
			if ($tokens[2] eq '-') {
				($from, $to) = ($to, $from);
			}
			$bsml_generator->AddExon($1, $from, $to);
		}
	}
	$bsml_generator->WriteBsml($out, $sequence_id, $project,
				   $gene_finder_name, $fasta);
}
