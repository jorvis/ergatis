#!/usr/local/bin/perl
use strict;
use warnings;

use Getopt::Long qw(:config no_ignore_case bundling);
use File::Basename;
use IO::File;

BEGIN {
use GenePredictionUtils::GeneMarkBsmlGenerator;
}

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

&parse_options;
&create_bsml;

sub print_usage
{
	my $progname = basename($0);
	die << "END";
usage: $progname [-input|i <input_genscan_data>] [-output|o <output_bsml>]
	[-fasta_file|f <fasta_data>] [-project|p <project_name>]
	[-help|h]
END
}

sub parse_options
{
	my %opts = ();
	GetOptions(\%opts, "input|i=s", "output|o=s", "fasta_file|f=s",
		   "project|p=s", "help|h");
	&print_usage if $opts{help};
	$in->open($opts{input}, "r") or
		die "Error reading genscan input $opts{input}: $!"
		if $opts{input};
	$out = $opts{output} if $opts{output};
	$fasta = $opts{fasta_file} if $opts{fasta_file};
	$project = $opts{project} if $opts{project};
}

sub create_bsml
{
	my $bsml_generator = new GenePredictionUtils::GeneMarkBsmlGenerator;
	my %gene_data = ();
	my $seq_id = undef;
	while (my $line = <$in>) {
		chomp $line;
		if ($line =~ /^Sequence name: ([\w\d\._]+)\s*/) {
			$seq_id = $1;
			next;
		}
		next if !length($line);
		$line =~ s/^\s+//g;
		my @tokens = split /\s+/, $line;
		next if $tokens[0] !~ /^(\d+)/;
		if (exists EXON_LABELS->{$tokens[3]}) {
			my $from = $tokens[4] - 1;
			my $to = $tokens[5] - 1;
			if ($tokens[2] eq '-') {
				($from, $to) = ($to, $from);
			}
			$bsml_generator->AddExon($1, $from, $to);
		}
	}
	$bsml_generator->WriteBsml($out, $seq_id, $project,
				   $gene_finder_name, $fasta);
}
