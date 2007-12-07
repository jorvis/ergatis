#!/usr/bin/perl

use strict;
use warnings;

use Getopt::Long qw(:config no_ignore_case bundling);
use XML::Twig;
use File::Basename;
use IO::File;

my $in	= undef;
my $out	= *STDOUT;

&parse_opts;

my %contigs = ();
my $twig = new XML::Twig(twig_roots =>
		{ 'Sequence[@class="contig"]' => \&sequence_handler } );
$twig->parsefile($in);

&print_agp;

sub parse_opts
{
	my %opts = ();
	GetOptions(\%opts, "input|i=s", "output|o=s", "help|h");
	&print_usage if $opts{help};
	$in = $opts{input} if $opts{input};
	$out = new IO::File($opts{output}, "w") or
		die "Error writing to output $opts{output}: $!"
		if $opts{output};
}

sub print_usage
{
	my $progname = basename($0);
	die << "END";
usage: $progname
END
}

sub sequence_handler
{
	my ($twig, $seq) = @_;
	if ($seq->class eq "contig") {
		my $contig_id = $seq->id;
		my $contig_len = $seq->att("length");
		my $numbering = $seq->first_child("Numbering");
		my ($from, $to, $strand);
		my $supercontig_id = $numbering->att("seqref");
		if ($numbering->att("ascending") == 1) {
			$from = $numbering->att("refnum") + 1;
			$to = $numbering->att("refnum") + $contig_len;
			$strand = "+";
		}
		else {
			$to = $numbering->att("refnum");
			$from = $numbering->att("refnum") - $contig_len + 1;
			$strand = "-";
		}
		push @{$contigs{$supercontig_id}},
			[$contig_id, $from, $to, $strand];
	}
	$twig->purge;
}

sub print_agp
{
	while (my ($supercontig_id, $contig_data) = each %contigs) {
		my @sorted_contig_data = sort { $a->[1] <=> $b->[1] }
						@{$contig_data};
		my $last = 1;
		my $part_num = 0;
		foreach my $contig (@sorted_contig_data) {
			my $id = $contig->[0];
			my $from = $contig->[1];
			my $to = $contig->[2];
			my $strand = $contig->[3];
			if ($from - 1 > $last) {
				print_gap($supercontig_id, $last, $from - 1, ++$part_num);
			}
			print_contig($supercontig_id, $contig, ++$part_num);
			$last = $to;
		}
	}
}

sub print_contig
{
	my ($supercontig_id, $contig, $part_num) = @_;
	my $id = $contig->[0];
	my $from = $contig->[1];
	my $to = $contig->[2];
	my $strand = $contig->[3];
	my $type = "W";
	my $contig_from = 1;
	my $contig_to = $to - $from + 1;
	$out->print("$supercontig_id\t$from\t$to\t$part_num\t$id\t$contig_from\t$contig_to\t$strand\n");
}

sub print_gap
{
	my ($supercontig_id, $from, $to, $part_num) = @_;
	my $type = "N";
	my $length = $to - $from;
	my $gap_type = "fragment";
	my $linkage = "yes";
	$out->print("$supercontig_id\t$from\t$to\t$part_num\t$type\t$length\t$gap_type\t$linkage\n");
}
