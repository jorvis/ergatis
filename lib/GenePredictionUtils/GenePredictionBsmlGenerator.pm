package GenePredictionUtils::GenePredictionBsmlGenerator;

use strict;
use warnings;
use BSML::BsmlBuilder;

## keep things properly scoped
{

my %_atts = ( 
    id_generator => undef,
);

sub new
{
	my ($class, %args) = @_;
    
	my $this = bless { %_atts }, $class;
	$this->Init;
    
    ## an id_generator reference is required
    if ( ! $args{id_generator} ) {
        die "you must pass an id_generator object reference when creating a new GenePredictionBsmlGenerator";
    }
    
    $this->{id_generator} = $args{id_generator};
    
	return $this;
}

sub Init
{
	my $this = shift;
	$this->{EXONS} = ();
}

sub AddExon
{
	my $this = shift;
	my ($gene_id, $from, $to) = @_;
	push @{$this->{EXONS}{$gene_id}}, [$from, $to];
}

sub WriteBsml
{
	my $this = shift;
	my $out = $_[0];
	my $doc = $this->GenerateBsml(@_);
	$doc->write($out);
}

sub GenerateBsml
{
	my $this = shift;
	my ($out, $seq_id, $project, $gene_finder_name, $fasta,
	    $feat_group_map) = @_;

	my $doc = new BSML::BsmlBuilder;
	my $idcreator = $this->{id_generator};
	my $seq = $doc->createAndAddSequence($seq_id, undef, '', 'dna', 'assembly');
	AddLink($seq, $gene_finder_name, 'input_of');
	my $feat_table = $doc->createAndAddFeatureTable($seq);
	if ($fasta) {
		$doc->createAndAddSeqDataImport($seq, 'fasta', $fasta,
						'', $seq_id);
	}
	while (my ($gene_id, $exons) = each %{$this->{EXONS}}) {
		my $gene_from = undef;
		my $gene_to = undef;
		my $gene_plus = 1;
		foreach my $coords (@$exons) {
			my $coord1 = $$coords[0];
			my $coord2 = $$coords[1];

			if ($coord1 > $coord2) {
				$gene_plus = 0;
				($coord1, $coord2) = ($coord2, $coord1);
			}
			if (!defined $gene_from || $coord1 < $gene_from) {
				$gene_from = $coord1;
			}
			if (!defined $gene_to || $coord2 > $gene_to) {
				$gene_to = $coord2;
			}
		}
		if (!$gene_plus) {
			($gene_from, $gene_to) = ($gene_to, $gene_from);
		}

		my $gene = $doc->createAndAddFeature
			($feat_table,
			 $idcreator->next_id(project		=> $project,
					    type	=> 'gene'),
			 '', 'gene');
             
		AddLink($gene, $gene_finder_name, 'computed_by');
        AddLoc($gene, $gene_from, $gene_to);

		my $feat_group = $doc->createAndAddFeatureGroup($seq, '', $gene->returnattr('id'));
		if ($feat_group_map) {
			$$feat_group_map{$gene_id} = $feat_group;
		}

		$feat_group->addBsmlFeatureGroupMember($gene->returnattr('id'), $gene->returnattr('class'));
        
        AddFeature($doc, $feat_table, $feat_group,
				   $idcreator, $project, 'CDS',
				   $gene_finder_name,
				   $gene_from, $gene_to);
        
        AddFeature($doc, $feat_table, $feat_group,
				   $idcreator, $project, 'transcript',
				   $gene_finder_name,
				   $gene_from, $gene_to);
        
		foreach my $coords (@$exons) {
			AddFeature($doc, $feat_table, $feat_group,
				   $idcreator, $project, 'exon',
				   $gene_finder_name,
				   $$coords[0], $$coords[1]);
		}
	}
	$doc->createAndAddAnalysis
		(id		=> "$gene_finder_name\_analysis",
		 sourcename	=> "$out",
         program => $gene_finder_name,
         version => 'current',
         algorithm => $gene_finder_name);
	return $doc;
}

sub AddLoc
{
	my ($feat, $from, $to) = @_;
	if ($from <= $to) {
		$feat->addBsmlIntervalLoc($from, $to, 0);
	}
	else {
		$feat->addBsmlIntervalLoc($to, $from, 1);
	}
}

sub AddLink
{
	my ($elt, $gene_finder_name, $role) = @_;
	$elt->addBsmlLink('analysis', "#$gene_finder_name\_analysis", $role);
}

sub AddFeature
{
	my ($doc, $feat_table, $feat_group, $idcreator, $project, $type,
	    $gene_finder_name, $coord1, $coord2) = @_;
	my $feat = $doc->createAndAddFeature
		($feat_table,
		 $idcreator->next_id(project		=> $project,
				    type	=> $type),
				 '', $type);
	AddLink($feat, $gene_finder_name, 'computed_by');
	$feat_group->addBsmlFeatureGroupMember($feat->returnattr('id'),
				 	       $feat->returnattr('class'))
		if $feat_group;
	AddLoc($feat, $coord1, $coord2);
}

}

1;
