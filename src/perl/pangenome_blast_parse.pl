#!/usr/local/bin/perl

use lib (@INC,$ENV{"PERL_MOD_DIR"});
no lib "$ENV{PERL_MOD_DIR}/i686-linux";
no lib ".";

=head1  NAME

pangenome_blast_parse.pl - Parses BLAST BSML output files and extracts data needed for pangenome analysis.

=head1 SYNOPSIS

USAGE: pangenome_query_list.pl
        --input=/path/to/somefile.bsml
		--output=/path/to/output.dat
      [ --log=/path/to/some/log ]

=head1 OPTIONS

B<--input_bsml,-i>
    BSML file containing BLAST results to parse. 

B<--output_path,-o>
    Desired prefix of output files.

B<--log,-d>
    optional. Will create a log file with summaries of all actions performed.

B<--debug>
    optional. Will display verbose status messages to STDERR if logging is disabled.

B<--help,-h>
    This help message/documentation.

=head1   DESCRIPTION

	This script parses the BSML output from WU-BLASTP and WU-TBLASTN runs for pangenome analysis
	and serializes the prepared data so that it can be later unserialized for the analysis step.

=head1 INPUT

    The input should be one of the set of BSML files containing BLAST results from a pangenome
	pipeline run.

=head1 OUTPUT

    There is no output unless you use the --log option.

=head1 CONTACT

    Brett Whitty
    bwhitty@tigr.org

=cut

use File::Basename;
use Pod::Usage;
use Storable;
use XML::Twig;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use strict;

my %options = ();
my $results = GetOptions (	\%options,
                          	'input|i=s',
              				'output_path|o=s',
#							'filter|f=s',
              				'debug|d=s',
              				'command_id=s',       ## passed by workflow
              				'logconf=s',          ## passed by workflow (not used)
              				'log|l=s',
              				'help|h',
		  				 ) || pod2usage();

my @results = ();
my @dups_temp = ();
my %dups = ();
my $filter;
my $skip_filter = 1;

my $twig = XML::Twig->new(
							twig_roots  => { 
												'Seq-pair-alignment' => \&processSeqPairAlignment
										   }
                         );

if (!-d $options{'output_path'}) {
	die "must specify output path as a directory";
}

$options{'output_path'} =~ s/\///;

if (-e $options{'output_path'}."/pangenome.filter.stored") {
	print STDERR "unserializing filter\n";
	$filter = retrieve($options{'output_path'}."/pangenome.filter.stored") || die "failed unserializing seq filter";
	$skip_filter = 0;
}

if (-e $options{'input'}) {
	$twig->parsefile($options{'input'});
} else {
	die "specified input file does not exist";
}

my $input_prefix = basename($options{'input'},".bsml");
#print $input_prefix."\n";

@dups_temp = sort(@dups_temp);
if (scalar(@dups_temp) > 1) {
	$dups{join(" ", @dups_temp)} = 1;
}

#store(\@results, $output_path."/".$input_prefix.".blast.stored") || die "couldn't serialize blast results array";
#store(\%dups, $output_path."/".$input_prefix.".dups.stored") || die "couldn't serialize dups hash";
store([\@results,\%dups], $options{'output_path'}."/".$input_prefix.".blast.stored") || die "couldn't serialize results";

exit(0);

sub processSeqPairAlignment { 
    my ($twig, $feat) = @_;
    my ($p_sim, $p_ident, $p_cov_ref, $p_cov_comp, $this_chain);
    
    my $query_id   = $feat->{'att'}->{'refseq'};
    my $subject_id = $feat->{'att'}->{'compseq'};

	my $len_total = 0;
	my $pident_sum = 0;
	my $psim_sum = 0;
	
    my $parent_chain = 0;
    
    for my $seqpairrun ( $feat->children('Seq-pair-run') ) {
		my $run_length = $seqpairrun->{att}->{runlength};
		$len_total += $run_length;
        for my $attribute ( $seqpairrun->children('Attribute') ) {
            if ($attribute->{att}->{'name'} eq 'percent_identity') {
                $p_ident = $attribute->{att}->{'content'};
            } elsif ($attribute->{att}->{'name'} eq 'percent_similarity') {
                $p_sim = $attribute->{att}->{'content'};
            } elsif ($attribute->{att}->{'name'} eq 'percent_coverage_refseq') {
                $p_cov_ref = $attribute->{att}->{'content'};
            } elsif ($attribute->{att}->{'name'} eq 'percent_coverage_compseq') {
                $p_cov_comp = $attribute->{att}->{'content'};
			}
        }
		$pident_sum += (($p_ident/100.0) * $run_length); 
		$psim_sum += (($p_sim/100.0) * $run_length); 
    }
	my $all_seg_p_ident = sprintf("%.1f", $pident_sum / $len_total * 100);
	my $all_seg_p_sim = sprintf("%.1f", $psim_sum / $len_total * 100);

	$query_id =~ /^([^_]+)_(.*)_[^_]+$/ || die "couldn't parse out db and sequence id from query\n$query_id";
	my($qdb, $qprot) = ($1, $2);
	$subject_id =~ /^([^_]+)_(.*)_[^_]+$/ || die "couldn't parse out db and sequence id from subject";
	my($sdb, $sprot) = ($1, $2);
	
	if ($skip_filter || ($filter->{$qdb}->{$qprot} && $filter->{$sdb}->{$sprot})) {
	#	if ($all_seg_p_ident >= 50.0 && ($p_cov_ref >= 50.0 || $p_cov_comp >= 50.0)) {  ## switched from p_ident to p_sim
	#		push (@results, [$qdb,$qprot,$sdb,$sprot,$all_seg_p_ident,$p_cov_ref]);
		if ($all_seg_p_sim >= 50.0 && ($p_cov_ref >= 50.0 || $p_cov_comp >= 50.0)) {
			push (@results, [$qdb,$qprot,$sdb,$sprot,$all_seg_p_sim,$p_cov_ref]);
		}
	#	if ($qdb eq $sdb && $qprot ne $sprot && $all_seg_p_indent == 100 && $p_cov_ref == 100 && $p_cov_comp == 100.0) {
		# 	I won't exclude the self hit here because we do want the query id in the group of dups
		if ($qdb eq $sdb && $all_seg_p_ident == 100 && $p_cov_ref == 100 && $p_cov_comp == 100) {
			push(@dups_temp, $sprot);
			print $sprot."\n";
			print $p_cov_ref." ".$p_cov_comp."\n";
		}
	}
}
