#!/usr/local/bin/perl
use lib (@INC,$ENV{"PERL_MOD_DIR"});
no lib "$ENV{PERL_MOD_DIR}/i686-linux";
no lib ".";
use strict;
use warnings;

use Getopt::Long qw(:config no_ignore_case bundling);
use File::Basename;
use IO::File;

BEGIN {
use GenePredictionUtils::GenscanBsmlGenerator;
}

use constant
{
	EXON_LABELS	=>	{	"Init"	=>	1,
					"Intr"	=>	1,
					"Term"	=>	1,
					"Sngl"	=>	1
				},
	PROMOTER_LABELS	=>	{	"Prom"	=>	1
				},
	POLYA_LABELS	=>	{	"PlyA"	=>	1
				}
};

my $in			= new IO::File->fdopen(fileno(STDIN), "r");
my $out			= "/dev/stdout";
my $fasta		= undef;
my $project		= "unknown";
my $gene_finder_name	=  "genscan";

&parse_options;
&create_bsml;

sub print_usage
{
	my $progname = basename($0);
	die << "END";
usage: $progname [-input|i <input_genscan_data>] [-output|o <output_bsml>]
	[-fasta_file|f <fasta_data>] [-project|p <project_name>]
	[-genscan_plus|g]
	[-help|h]
END
}

sub parse_options
{
	my %opts = ();
	GetOptions(\%opts, "input|i=s", "output|o=s", "fasta_file|f=s",
		   "project|p=s", "genscan_plus|g", "help|h");
	&print_usage if $opts{help};
	$gene_finder_name = "genscan_plus" if $opts{genscan_plus};
	$in->open($opts{input}, "r") or
		die "Error reading genscan input $opts{input}: $!"
		if $opts{input};
#	$out->open($opts{output}, "w") or
#		die "Error writing BSML output $opts{output}: $!"
#		if $opts{output};
	$out = $opts{output} if $opts{output};
	$fasta = $opts{fasta_file} if $opts{fasta_file};
	$project = $opts{project} if $opts{project};
}

sub create_bsml
{
	my $bsml_generator = new GenePredictionUtils::GenscanBsmlGenerator;
	my %gene_data = ();
	my $seq_id = undef;
	while (my $line = <$in>) {
		chomp $line;
		if ($line =~ /^Sequence ([\w\d\._]+)\s/) {
			$seq_id = $1;
			next;
		}
		next if !length($line);
		$line =~ s/^\s+//g;
		my @tokens = split /\s+/, $line;
		next if $tokens[0] !~ /^(\d+)\.\d+/;

		my $from = $tokens[3] - 1;
		my $to = $tokens[4] - 1;

		if (exists EXON_LABELS->{$tokens[1]}) {
			$bsml_generator->AddExon($1, $from, $to);
		}
		elsif (exists PROMOTER_LABELS->{$tokens[1]}) {
			$bsml_generator->AddPromoter($1, $from, $to);
		}
		elsif (exists POLYA_LABELS->{$tokens[1]}) {
			$bsml_generator->AddPolyA($1, $from, $to);
		}
	}
	$bsml_generator->WriteBsml($out, $seq_id, $project,
				   $gene_finder_name, $fasta);
}
