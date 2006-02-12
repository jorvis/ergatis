#!/usr/local/bin/perl
use strict;
use warnings;

use Getopt::Long qw(:config no_ignore_case bundling);
use File::Basename;
use IO::File;

BEGIN {
use GenePredictionUtils::FGeneShBsmlGenerator;
}

use constant
{
	EXON_LABELS	=>	{	"CDSf"	=>	1,
					"CDSi"	=>	1,
					"CDSl"	=>	1,
					"CDSo"	=>	1
				},
	TSS_LABELS	=>	{	"TSS"	=>	1
				},
	POLYA_LABELS	=>	{	"PolA"	=>	1
				}
};

my $in			= new IO::File->fdopen(fileno(STDIN), "r");
my $out			= "/dev/stdout";
my $fasta		= undef;
my $project		= "unknown";
my $gene_finder_name	= "fgenesh";

&parse_options;
&create_bsml;

sub print_usage
{
	my $progname = basename($0);
	die << "END";
usage: $progname [-input|i <input_fgenesh_data>] [-output|o <output_bsml>]
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
	my $bsml_generator = new GenePredictionUtils::FGeneShBsmlGenerator;
	my %gene_data = ();
	my $seq_id = undef;
	while (my $line = <$in>) {
		chomp $line;
		$line =~ s/^\s+//g;
		if ($line =~ /^Seq name: ([\w\d\._]+)\s/) {
			$seq_id = $1;
			next;
		}
		next if !length($line);
		my @tokens = split /\s+/, $line;
		next if $tokens[0] !~ /(\d+)/;
		if (scalar(@tokens) == 12 &&
		    exists EXON_LABELS->{$tokens[3]}) {
			my $from = $tokens[4] - 1;
			my $to = $tokens[6] - 1;
			($from, $to) = ($to, $from) if
				$tokens[1] eq '-';
			$bsml_generator->AddExon($1, $from, $to);
		}
		elsif (scalar(@tokens) == 5) {
			my $from = $tokens[3] - 1;
			my $to = $from + 1;
			($from, $to) = ($to, $from) if
				$tokens[1] eq '-';
			if (exists TSS_LABELS->{$tokens[2]}) {
				$bsml_generator->AddTSS($1, $from, $to);
			}
			elsif (exists POLYA_LABELS->{$tokens[2]}) {
				$bsml_generator->AddPolyA($1, $from, $to);
			}
		}
	}
	$bsml_generator->WriteBsml($out, $seq_id, $project,
				   $gene_finder_name, $fasta);
}
