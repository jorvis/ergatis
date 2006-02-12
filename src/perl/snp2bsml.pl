#!/usr/local/bin/perl
use lib (@INC,$ENV{"PERL_MOD_DIR"});
no lib "$ENV{PERL_MOD_DIR}/i686-linux";
no lib ".";
=head1 NAME

snp2bsml.pl - convert clustered snp output information from show-snps to bsml

=head1 SYNOPSIS

usage: snp2bsml.pl [ --input|-i=/path/to/show-snps_clustered_output ]
		   [ --output|-o=/path/to/snp_data.bsml ]
		     --database|-d=database_to_be_loaded
		   [ --log|-l=/path/to/logfile ]
		   [ --debug|-d=debug_level ]
		   [ --moltype|-m=molecule_type ]
		   [ --seqdata|-s=sequence_data ]

=head1 OPTIONS

B<--input, -i>
	Clustered output from show-snps tool from MUMer suite

B<--output, -o>
	BSML multiple sequence alignment file for SNP representation

B<--database, -d>
	Database to load BSML data into

B<--log, -l>
	Log file

B<--debug, -D>
	Debug level

B<--moltype, -m>
	Molecule type (nucleotide, protein)
	Defaults to nucleotide

B<--seqdata, -s>
	SNP sequence data, in tab delimited format, containing:
	sequence_id, position, orientation, base_call

B<--help, -h>
	This help message

=head1 DESCRIPTION

This script will parse out clustered output from show-snps tool from the
MUMer package and create a multiple sequence alignment BSML representation
for the SNP data.

=head1 INPUT

The input for this script is the output from show-snps tool from the MUMer
package, in tab delimited format.

=head1 OUTPUT

The output for this script is the BSML multiple sequence alignment
representation of the SNP data.

=head1 CONTACT

Ed Lee (elee@tigr.org)

=cut

use strict;
use warnings qw(all);

use IO::File;
use BSML::BsmlDoc;
use Papyrus::TempIdCreator;
use Workflow::Logger;
use Pod::Usage;
use Getopt::Long qw(:config no_ignore_case);
use Sequence::SeqLoc;

my $moltype		= 'nucleotide';
my $in			= undef;
my $out			= undef;
my $db			= undef;
my $logger		= undef;
my %seqdata		= ();
my $RECORD_DELIM	= '//';

&get_opts;
&parse_data;

sub print_usage
{
	pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDERR} );
}

sub get_opts
{
	my %opts	= ();
	my $input_file	= '/dev/stdin';
	my $output_file	= '/dev/stdout';
	my $log_file	= Workflow::Logger::get_default_logfilename;
	my $debug	= 4;
	GetOptions(\%opts,
		   'in|i=s', 'out|o=s', 'database|d=s',
		   'log|l=s', 'debug|D=i', 'moltype|m=s', 'seqdata|s=s',
		   'help|h');
	print_usage if $opts{help};
	$input_file = $opts{in} if $opts{in};
	$output_file = $opts{out} if $opts{out};
	$db = $opts{database} if $opts{database};
	$log_file = $opts{log} if $opts{log};
	$debug = $opts{debug} if $opts{debug};
	$logger = new Workflow::Logger('LOG_FILE'	=> $log_file,
				       'LOG_LEVEL'	=> $debug)->
		  get_logger;
	$in = new IO::File($input_file, "r") or
		$logger->logdie("Error reading input $input_file: $!");
	$out = new IO::File($output_file, "w") or
		$logger->logdie("Error writing output $output_file: $!");
	$logger->logdie("No database provided") if !$db;
	$moltype = $opts{moltype} if $opts{moltype};
	init_seqdata($opts{seqdata}) if $opts{seqdata};
	$logger->logdie("No SNP seqdata provided")
		if (!scalar(keys(%seqdata)));
}

sub init_seqdata
{
	my $seqdata_fname = shift;
	my $data = new IO::File($seqdata_fname) or
		die "Error reading SNP sequence data $seqdata_fname: $!";
	while (my $line = <$data>) {
		chomp $line;
		my @tokens = split /\t/, $line;
		$seqdata{$tokens[0]}{$tokens[1]}{$tokens[2]} = $tokens[3];
	}
}

sub parse_data
{
	my $id_creator		= new Papyrus::TempIdCreator;
	my $doc			= new BSML::BsmlDoc;
	my @indels		= ();
	my @cluster		= ();
	while (my $line = <$in>) {
		chomp $line;
		if ($line =~ /^$RECORD_DELIM$/) {
			add_table($doc, \@cluster, $id_creator, 'SNP');
			if (is_indel(\@cluster)) {
				push @indels, [];
				push @{$indels[$#indels]}, @cluster;
			}
			else {
				add_table($doc, \@cluster, $id_creator,
					  'indel');
			}
			@cluster = ();
			next;
		}
		my $snp = new Sequence::SeqLoc($line);
		$snp->SetSeqData(get_seq_data($snp)) if $snp->GetLength; 
		push @cluster, $snp;
	}
	process_indels($doc, \@indels, $id_creator);
	$doc->write($out);
}

sub add_table
{
	my ($doc, $cluster, $id_creator, $type) = @_;
	my $msa_table = $doc->returnBsmlMultipleAlignmentTableR
			($doc->addBsmlMultipleAlignmentTable);
	$msa_table->addBsmlLink('analysis', '#snp_analysis', 'computed_by');
	$msa_table->addattr('molecule-type', $moltype);
	$msa_table->addattr('id', $id_creator->new_id(db => $db,
					   	      so_type => $type));
	$msa_table->addattr('class', $type);
	my $aln_summary = $msa_table->returnBsmlAlignmentSummaryR
			  ($msa_table->addBsmlAlignmentSummary);
	$aln_summary->addattr('seq-format', 'msa');
	$aln_summary->addattr('seq-type', $moltype);
	my $seq_align = $msa_table->returnBsmlSequenceAlignmentR
			($msa_table->addBsmlSequenceAlignment);
	my $seq_nums = "";
	for (my $i = 0; $i < scalar(@{$cluster}); ++$i) {
		$seq_nums .= "$i";
		$seq_nums .= ":" if $i < $#$cluster;
		my $loc = $$cluster[$i];
		my $aln_seq = $aln_summary->returnBsmlAlignedSequenceR
			      ($aln_summary->addBsmlAlignedSequence);
		$aln_seq->addattr('seqref', $loc->GetId);
		$aln_seq->addattr('length',
			$loc->GetSeqData ? length($loc->GetSeqData) :
			0);
		$aln_seq->addattr('seqnum', $i);
		$aln_seq->addattr('start', $loc->GetFrom);
		$aln_seq->addattr('on-complement',
				  $loc->GetStrand eq '+' ? '0' : '1');
		my $seq_data = $seq_align->returnBsmlSequenceDataR
			       ($seq_align->addBsmlSequenceData);
		$seq_data->addattr('seq-name', $loc->GetId);
		$seq_data->addSequenceAlignmentData($loc->GetSeqData)
			if $loc->IsSetSeqData;
	}
	$seq_align->addattr('sequences', $seq_nums);
}

sub process_indels
{
	my ($doc, $indels, $id_creator) = @_;
	my @sorted_clusters = sort cluster_comparator @{$indels};
	my @adjacent_clusters = ();
	my %written_clusters = ();
	my $prev_adjacent = 0;
	for (my $i = 0; $i < scalar(@sorted_clusters) - 1; ++$i) {
		my $cluster1 = $sorted_clusters[$i];
		for (my $j = $i + 1; $j < scalar(@sorted_clusters); ++$j) {
			my $cluster2 = $sorted_clusters[$j];
			if (is_adjacent($cluster1, $cluster2)) {
				if (!$prev_adjacent) {
					push @adjacent_clusters, [];
					push @{$adjacent_clusters
						[$#adjacent_clusters]},
						$cluster1;
					++$written_clusters{$cluster1};
					$prev_adjacent = 1;
				}
				push @{$adjacent_clusters
					[$#adjacent_clusters]},
					$cluster2;
				++$written_clusters{$cluster2};
				last;
			}
			else {
				$prev_adjacent = 0;
			}
		}
	}
	foreach my $cluster_group (@adjacent_clusters) {
		my $indel = "";
		my $loc_with_seq = undef;
		foreach my $cluster (@{$cluster_group}) {
			foreach my $loc (@{$cluster}) {
				if ($loc->IsSetSeqData) {
					$indel .= $loc->GetSeqData;
					$loc_with_seq = $loc
						if !$loc_with_seq;
				}
			}
		}
		my @indel_cluster = ();
		foreach my $loc (@{$$cluster_group[0]}) {
			if ($loc->GetId eq $loc_with_seq->GetId) {
				my $indel_loc =
					create_new_loc($loc, $indel);
				push @indel_cluster, $indel_loc;
			}
			else {
				push @indel_cluster, $loc;
			}
		}
		add_table($doc, \@indel_cluster, $id_creator, 'indel');
	}
	foreach my $cluster (@sorted_clusters) {
		next if exists $written_clusters{$cluster};
		add_table($doc, \@{$cluster}, $id_creator, 'indel');
	}
}

sub copy_snp_data
{
	my ($src, $dest) = @_;
	$dest->SetQueryId($src->GetQueryId);
	$dest->SetQueryPosition($src->GetQueryPosition);
	$dest->SetQueryFrame($src->GetQueryFrame);
	$dest->SetSubjectId($src->GetSubjectId);
	$dest->SetSubjectPosition($src->GetSubjectPosition);
	$dest->SetSubjectFrame($src->GetQueryFrame);
}

sub snp_comparator
{
	my ($snp1, $snp2, $is_query) = @_;
	if ($is_query) {
		if ($snp1->GetQueryFrame > 0) {
			return $snp1->GetQueryPosition <
				$snp2->GetQueryPosition;
		}
		else {
			return $snp1->GetQueryPosition >
				$snp2->GetQueryPosition;
		}
	}
	else {
		if ($snp1->GetSubjectFrame > 0) {
			return $snp1->GetSubjectPosition <
				$snp2->GetSubjectPosition;
		}
		else {
			return $snp1->GetSubjectPosition >
				$snp2->GetSubjectPosition;
		}
	}
}

sub is_indel
{
	my $cluster = shift;
	foreach my $loc (@{$cluster}) {
		return 1 if !$loc->GetLength;
	}
	return 0;
}

sub cluster_comparator
{
	if ($$a[0]->GetId ne $$b[0]->GetId) {
		return $$a[0]->GetId cmp $$b[0]->GetId;
	}
	my $min_cluster_size = scalar(@{$a}) < scalar(@{$b}) ?
		scalar(@{$a}) : scalar(@{$b});
	for (my $i = 0; $i < $min_cluster_size; ++$i) {
		if ($$a[$i]->GetId ne $$b[$i]->GetId) {
			return $$a[$i]->GetId cmp $$b[$i]->GetId;
		}
		if ($$a[$i]->GetFrom != $$b[$i]->GetFrom) {	
			return $$a[0]->GetFrom <=> $$b[0]->GetFrom;
		}
	}
}

sub is_adjacent
{
	my ($cluster1, $cluster2) = @_;
	return 0 if scalar(@{$cluster1}) != scalar(@{$cluster2});
	for (my $i = 0; $i < scalar(@{$cluster1}); ++$i) {
		my $loc1 = $$cluster1[$i];
		my $loc2 = $$cluster2[$i];
		if ($loc1->GetId ne $loc2->GetId) {
			return 0;
		}
		if ($loc1->GetFrom != $loc2->GetFrom &&
		    $loc1->GetFrom + 1 != $loc2->GetFrom) {
			return 0;
		}
	}
	return 1;
}

sub get_seq_data
{
	my $loc = shift;
	return $seqdata{$loc->GetId}{$loc->GetFrom}{$loc->GetStrand};
}

sub create_new_loc
{
	my ($loc, $indel) = @_;
	my $new_loc = new Sequence::SeqLoc($loc);
	$new_loc->SetSeqData($indel);
#	if ($new_loc->GetStrand eq '+') {
#		$new_loc->SetTo($new_loc->GetFrom + length($indel));
#	}
#	else {
#		$new_loc->SetFrom($new_loc->GetTo - length($indel));
#	}
	return $new_loc;
}
