=head1 NAME

ber_blast_hit_selector.pl - flter hits from blast btab output

=head1 SYNOPSIS

usage: ber_blast_hit_selector.pl

	[
		--input|i=/path/to/blast_btab_data
		--output|o=/path/to/filtered_blast_btab_data
		--min_bit_score|b=min_bit_score
		--min_raw_score|r=min_raw_score
		--min_pct_id|p=min_percent_identity
		--min_pct_sim|P=min_percent_similarity
		--max_eval|e=max_e_value
		--max_pval|E=max_p_value
		--max_num_hits|n=max_number_of_hits_to_output
		--max_num_hits_per_region|N=max_number_of_hits_per_region
		--min_num_experimental|V=min_number_experimental
		--help|h
	]

=head1 OPTIONS

B<--input, i>
	BLAST btab data [default = stdin]

B<--output, o>
	Filtered BLAST btab data [default = stdout]

B<--min_bit_score, b>
	Minimum bitscore for hit [default = 0]

B<--min_raw_score, r>
	Minimum rawscore for hit [default = 0]

B<--min_pct_id, p>
	Minimum percent identity for hit [default = 0]

B<--min_pct_sim, P>
	Minimum percent similarity for hit [default = 0]

B<--max_eval, e>
	Maximum e-value for hit [default = 1]

B<--max_pval, E>
	Maximum p-value for hit [default = 1]

B<--max_num_hits, n>
	Maximum number of hits to output [default = INT_MAX]

B<--max_num_hits_per_region, N>
	Maximum number of hits to output per region [default = INT_MAX]

B<--min_num_experimental, V>
	Minimum number of hits that are experimentally verified
	[default = 0]

B<--help, h>
	This help screen

=head1 DESCRIPTION

This script will filter out hits from BLAST btab output (given the parameters)
and display the n best hits on bitscore (where n is the max number of hits to
display).

=head1 CONTACT

Ed Lee (elee@tigr.org)

=cut

use strict;
use warnings;

use Getopt::Long qw(:config no_ignore_case bundling);
use IO::File;
use Pod::Usage;
use Blast::BlastHitDataType;
use POSIX;

my $min_bit_score		= 0;
my $min_raw_score		= 0;
my $min_pct_id			= 0;
my $min_pct_sim			= 0;
my $max_eval			= 1;
my $max_pval			= 1;
my $max_num_hits		= POSIX::INT_MAX;
my $max_num_hits_per_region	= 0;
my $min_num_experimental	= 0;
my $in				= new IO::File->fdopen(0, "r");
my $out				= new IO::File->fdopen(1, "w");
my @hits			= ();

&get_options;
&filter_hits;
&print_hits;

sub print_usage
{
	pod2usage( {-exitval => 1, -verbose => 2, -output => \*STDOUT} );
}

sub get_options
{
	my %opts;
	GetOptions(\%opts,
		   "input|i=s", "output|o=s",
		   "min_bit_score|b=f", "min_raw_score|r=i",
		   "min_pct_id|p=f", "min_pct_sim|P=f",
		   "max_eval|e=f", "max_pval|E=f",
		   "max_num_hits|n=i", "max_num_hits_per_region|N=i",
		   "min_num_experimental|V=i",
		   "help|h") or &print_usage;
	&print_usage if $opts{help};
	$in->open($opts{input}, "r") if $opts{input};
	$out->open($opts{output}, "w") if $opts{output};
	$min_bit_score = $opts{min_bit_score} if defined $opts{min_bit_score};
	$min_raw_score = $opts{min_raw_score} if defined $opts{min_raw_score};
	$min_pct_id = $opts{min_pct_id} if defined $opts{min_pct_id};
	$min_pct_sim = $opts{min_pct_sim} if defined $opts{min_pct_sim};
	$max_eval = $opts{max_eval} if defined $opts{max_eval};
	$max_pval = $opts{max_pval} if defined $opts{max_pval};
	$max_num_hits = $opts{max_num_hits} if defined $opts{max_num_hits};
	$max_num_hits_per_region = $opts{max_num_hits_per_region}
		if defined $opts{max_num_hits_per_region};
	$min_num_experimental = $opts{min_num_experimental}
		if defined $opts{min_num_experimental};
}

sub filter_hits
{
	while (my $line = <$in>) {
		chomp $line;
		my $hit = new Blast::BlastHitDataType($line);
		if ($hit->GetBitScore >= $min_bit_score &&
		    $hit->GetRawScore >= $min_raw_score &&
		    $hit->GetPercentId >= $min_pct_id &&
		    $hit->GetPercentSimilarity >= $min_pct_sim &&
		    $hit->GetEValue <= $max_eval &&
		    $hit->GetPValue <= $max_pval) {
			push @hits, $hit;
		}
	}
	my @sorted_hits = sort score_comparator @hits;
	@hits = @sorted_hits;
}

sub print_hits
{
	if (!$max_num_hits_per_region) {
		print_filtered_hits(\@hits, $max_num_hits);
	}
	else {
		my @all_hits_by_region = ();
		foreach my $hit (@hits) {
			add_hits_by_region(\@all_hits_by_region, $hit);
		}
		foreach my $hits_by_region (@all_hits_by_region) {
			print_filtered_hits($hits_by_region,
					    $max_num_hits_per_region);
		}
	}
}

sub score_comparator
{
	return $b->GetBitScore <=> $a->GetBitScore
		if $a->GetBitScore > 0 && $b->GetBitScore > 0;
	return $b->GetRawScore <=> $a->GetRawScore;
}

sub add_hits_by_region
{
	my ($all_hits_by_region, $new_hit) = @_;
	foreach my $hits_by_region (@$all_hits_by_region) {
		foreach my $hit (@$hits_by_region) {
			if (overlaps($hit, $new_hit)) {
				push @$hits_by_region, $new_hit;
				return;
			}
		}

	}
	push @$all_hits_by_region, [$new_hit];
}

sub overlaps
{
	my ($hit1, $hit2) = @_;
	if ($hit1->GetQueryName ne $hit2->GetQueryName) {
		return 0;
	}
	if ($hit1->GetQueryBegin < $hit2->GetQueryBegin) {
		return $hit1->GetQueryEnd >= $hit2->GetQueryBegin;
	}
	else {
		return $hit2->GetQueryEnd >= $hit1->GetQueryBegin;
	}
}

sub print_filtered_hits
{
	my ($hits, $num_to_print) = @_;
	my $num_experimental = 0;
	for (my $i = 0; $i < scalar(@$hits); ++$i) {
		last if $i >= $num_to_print &&
			$num_experimental >= $min_num_experimental;
		my $ok = 1;
		if (is_experimentally_verified($$hits[$i])) {
			++$num_experimental;
		}
		elsif ($i >= $num_to_print) {
			$ok = 0;
		}
		print $out $$hits[$i]->ToString, "\n" if $ok;
	}
}

sub is_experimentally_verified
{
	my $hit = shift;
	return $hit->GetDescription =~ /experimental=1/;
}
