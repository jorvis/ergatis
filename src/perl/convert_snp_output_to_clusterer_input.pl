#!/usr/local/bin/perl
BEGIN{foreach (@INC) {s/\/usr\/local\/packages/\/local\/platform/}};
use lib (@INC,$ENV{"PERL_MOD_DIR"});
no lib "$ENV{PERL_MOD_DIR}/i686-linux";
no lib ".";
use strict;
use warnings;
use IO::File;
use Getopt::Long qw(:config no_ignore_case);

use MUMmer::SnpDataType;

my $in	= new IO::File->fdopen(fileno(STDIN), "r");
my $out	= new IO::File->fdopen(fileno(STDOUT), "w");

&parse_options;
&convert_output;

sub print_usage
{
	die << "END";
usage: convert_snp_output_to_clusterer_input.pl [-input|i <show-snps_output>]
	[-output|o <output>] [-help|h]
END
}

sub parse_options
{
	my %opts = ();
	GetOptions(\%opts, "input|i=s", "output|o=s", "help|h");
	print_usage if $opts{help};
	$in->open($opts{input}, "r") or
		die "Error reading from $opts{input}: $!"
		if $opts{input};
	$out->open($opts{output}, "w") or
		die "Error writing to $opts{output} : $!"
		if $opts{output};
}

sub convert_output
{
	while (my $line = <$in>) {
		chomp $line;
		if ($line =~ /^\[/) {
			MUMmer::SnpDataType::SetHeader($line);
			next;
		}
		my @tokens = split /\t/, $line;
		next if scalar(@tokens) !=
			MUMmer::SnpDataType::GetNumberOfColumns;
		my $snp = new MUMmer::SnpDataType($line);
		my $query_id = $snp->GetQueryId;
		my $subj_id = $snp->GetSubjectId;

		my $query_from;
		my $query_to;
		my $subj_from;
		my $subj_to;
		my $query_strand;
		my $subj_strand;
		&set_loc_data(\$query_from, \$query_to, \$query_strand,
			      $snp->GetQueryFrame,
			      $snp->GetQueryPosition,
			      $snp->GetQuerySubstitution);
		&set_loc_data(\$subj_from, \$subj_to, \$subj_strand,
			      $snp->GetSubjectFrame,
			      $snp->GetSubjectPosition,
			      $snp->GetSubjectSubstitution);
		&print_data($query_id, $query_from, $query_to, $query_strand);
		&print_data($subj_id, $subj_from, $subj_to, $subj_strand);
		$out->print("//\n");
	}
}

sub print_data
{
	my ($id, $from, $to, $strand) = @_;
	$out->printf("%s\t%d\t%d\t%s\n", $id, $from, $to, $strand);
}

sub set_loc_data
{
	my ($from, $to, $strand, $frame, $pos, $subst) = @_;
	#if ($frame > 0) {
		$$from = $pos;
		$$to = $subst  ne "." ?  $pos + 1 : $pos;
		$$strand = $frame > 0 ? "+" : "-";
	#}
	#else {
	#	$$from = $subst ne "." ?  $pos - 1 : $pos;
	#	$$to = $pos;
	#	$$strand = "-";
	#}
}
