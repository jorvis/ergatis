#!/usr/local/bin/perl
BEGIN{foreach (@INC) {s/\/usr\/local\/packages/\/local\/platform/}};
use lib (@INC,$ENV{"PERL_MOD_DIR"});
no lib "$ENV{PERL_MOD_DIR}/i686-linux";
no lib ".";

use strict;
use warnings;

use Getopt::Long qw(:config no_ignore_case bundling);
use IO::File;
use File::Basename;

use constant
{
	MASK_AMBIGUOUS	=> 0x1
};

my $in		= new IO::File->fdopen(fileno(STDIN), "r");
my $out		= new IO::File->fdopen(fileno(STDOUT), "w");
my $mask_flag	= 0;

&parse_options;
&mask_fasta;

sub print_usage
{
	my $progname = basename($0);
	die << "END";
usage: $progname [-input|i <fasta_input>] [-output|o <fasta_output>]
	[-mask_ambiguous|m]

	mask_ambiguous|m: mask ambiguous nucleotide bases
			  (non A C G T)
END
}

sub parse_options
{
	my %opts = ();
	GetOptions(\%opts, "input|i=s", "output|o=s",
		   "mask_ambiguous|m", "help|h") or &print_usage;
	&print_usage if $opts{help};
	$in->open($opts{input}, "r") 
		or die "Error reading fasta $opts{input}: $!"
		if $opts{input};
	$out->open($opts{output}, "w")
		or die "Error writing to $opts{output}: $!"
		if $opts{output};
	$mask_flag |= MASK_AMBIGUOUS if $opts{mask_ambiguous};
}

sub mask_fasta
{
	while (my $line = <$in>) {
		chomp $line;
		if ($line =~ /^>/) {
			print $out "$line\n";
		}
		else {
			&mask_sequence($line);
		}
	}
}

sub mask_sequence
{
	my $seq = uc(shift);
	if ($mask_flag & MASK_AMBIGUOUS) {
		#$seq =~ s/[^ACGT]/N/g;
		$seq =~ tr/ACGT/N/c;
	}
	print $out "$seq\n";
}
