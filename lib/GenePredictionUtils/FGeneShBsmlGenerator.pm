package GenePredictionUtils::FGeneShBsmlGenerator;
@ISA = qw(GenePredictionUtils::GenePredictionBsmlGenerator);

use strict;
use warnings;

BEGIN {
    require '/usr/local/devel/ANNOTATION/cas/lib/site_perl/5.8.5/GenePredictionUtils/GenePredictionBsmlGenerator.pm';
    import GenePredictionUtils::GenePredictionBsmlGenerator;
}

sub new
{
	my $class = shift;
	my $this = {};
	bless $this, $class;
	$this->Init;
	return $this;
}

sub Init
{
	my $this = shift;
	$this->{TSS} = ();
	$this->{POLYAS} = ();
	$this->SUPER::Init();
}

sub AddTSS
{
	my $this = shift;
	my ($gene_id, $from, $to) = @_;
	push @{$this->{TSS}{$gene_id}}, [$from, $to];
}

sub AddPolyA
{
	my $this = shift;
	my ($gene_id, $from, $to) = @_;
	push @{$this->{POLYAS}{$gene_id}}, [$from, $to];
}

sub GenerateBsml
{
	my $this = shift;
	my ($out, $seq_id, $project, $gene_finder_name, $fasta) = @_;
	my %feat_group_map = ();
	my $idcreator = new Papyrus::TempIdCreator;
	my $doc = $this->SUPER::GenerateBsml(@_, \%feat_group_map);
	my $seq = $doc->returnBsmlSequenceR(0);
	my $feat_table = $seq->returnBsmlFeatureTableR(0);
	my $add_feature_func =
		\&GenePredictionUtils::GenePredictionBsmlGenerator::AddFeature;
	while (my ($gene_id, $promoters) = each %{$this->{TSS}}) {
		foreach my $coords (@$promoters) {
			my $feat_group = $feat_group_map{$gene_id};
			$add_feature_func->
				($doc, $feat_table, $feat_group,
				 $idcreator, $project,
				 'transcription_start_site',
				 $gene_finder_name, $$coords[0], $$coords[1]);
		}
	}
	while (my ($gene_id, $polyas) = each %{$this->{POLYAS}}) {
		foreach my $coords (@$polyas) {
			my $feat_group = $feat_group_map{$gene_id} or next;
			$add_feature_func->
				($doc, $feat_table, $feat_group,
				 $idcreator, $project, 'polyA_signal_sequence',
				 $gene_finder_name, $$coords[0], $$coords[1]);
		}
	}
	return $doc;
}

1;
