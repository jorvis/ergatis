#!/usr/local/bin/perl
use lib (@INC,$ENV{"PERL_MOD_DIR"});
no lib "$ENV{PERL_MOD_DIR}/i686-linux";
no lib ".";
use strict;
use warnings;

use File::Basename;

use Getopt::Long qw(:config no_ignore_case);
use IO::File;
use MUMer::SnpDataType;

my $in		= new IO::File->fdopen(fileno(STDIN), "r");
my $out		= new IO::File("|sort|uniq>/dev/stdout");
my $INDEL_CHAR	= '.';

&parse_options;
&extract_data;

sub print_usage
{
	my $progname = basename($0);
	die << "END";
usage: $progname [--input|-i <show-snps_output>]
	[--output|-o <output>] [--help|-h]
END
}

sub parse_options
{
	my %opts = ();
	GetOptions(\%opts, "input|i=s", "output|o=s", "help|h");
	print_usage if $opts{help};
	$in->open($opts{input}, "r")
		or die "Error reading from $opts{input}: $!"
		if $opts{input};
	$out->open("|sort|uniq>$opts{output}")
		or die "Error writing to $opts{output}: $!"
		if $opts{output};
}

sub extract_data
{
	while (my $line = <$in>) {
		chomp $line;
		if ($line =~ /^\[/) {
			MUMer::SnpDataType::SetHeader($line);
			next;
		}
		my @tokens = split /\t/, $line;
		next if scalar(@tokens) !=
			MUMer::SnpDataType::GetNumberOfColumns;
		my $snp = new MUMer::SnpDataType($line);
		$out->printf("%s\t%d\t%s\t%s\n", $snp->GetQueryId,
			     $snp->GetQueryPosition,
			     $snp->GetQueryFrame > 0 ? "+" : "-",
			     $snp->GetQuerySubstitution)
			if $snp->GetQuerySubstitution ne $INDEL_CHAR;
		$out->printf("%s\t%d\t%s\t%s\n", $snp->GetSubjectId,
			     $snp->GetSubjectPosition,
			     $snp->GetSubjectFrame > 0 ? "+" : "-",
			     $snp->GetSubjectSubstitution)
			if $snp->GetSubjectSubstitution ne $INDEL_CHAR;
	}
}
