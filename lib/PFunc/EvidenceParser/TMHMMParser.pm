package PFunc::EvidenceParser::TMHMMParser;

use strict;
use warnings;
use XML::Twig;
use Data::Dumper;

use base qw(PFunc::EvidenceParser);

########## Class Vars ##########
my $annotation_type = "TMHMM";
my $transmembrane = {
    'gene_product_name' => "putative membrane protein",
    'gene_symbol'       => "",
    'EC'                => "",
    'GO'                => ['GO:0008150', 'GO:0003674', 'GO:0016020'],
    'TIGR_Role'         => 88,
};
my $source = "TMHMM";
my $minimum_membrane_spanning_regions = 3;
################################

sub new {
    my ($class, %args) = @_;
    my $self = $class->SUPER::new( %args );
    bless($self, $class);
    return $self;
}

sub _parse {
    my ($self, $fh) = @_;
    
    my $twig = new XML::Twig( 'twig_roots' => {
        'Sequence[@class="polypeptide"]' => sub { $self->_handle_polypeptide( @_ ) },
    } );
    $twig->parse( $fh );
}

sub _handle_polypeptide {
    my ($self, $twig, $el) = @_;
    
    my $annot_id = $self->lookup_feature_id( $el->att('id'), $self->annotate_on );
    my $annotation = $self->get_feature_annotation( $annot_id );
    return unless( $annotation );

    my $annot_flag = 0;

    #every file should have this. Even if no transmembrane regions exist.
    my ($tmh_count_att) = $el->find_nodes('Attribute[@name="tmh_count"]');
    my $tmh_count;
    if( $tmh_count_att ) {
        $tmh_count = $tmh_count_att->att('content');
    } else {
        die('Could not find Attribute[@name="tmh_count"] for sequence '.$el->att('id'));
    }

    if( $tmh_count >= $minimum_membrane_spanning_regions ) {
        $annot_flag = 1;
    } elsif( $tmh_count > 0 ) {
        
        #get the sequence length
        my $seq_length = ($el->find_nodes('Attribute[@name="length"]'))[0]->att('content');
        
        #get predicted aa's in transmembrane region
        my $exp_aa_in_tmh = ($el->find_nodes('Attribute[@name="exp_aa_in_tmh"]'))[0]->att('content');

        my $perc_in_tmh = int( ($exp_aa_in_tmh/$seq_length) * 10000 ) / 100;

        $annot_flag = 1 if( $perc_in_tmh > 50 );
    }

    if( $annot_flag ) {
        $self->_annotate_feature( $annotation );
    }
      
}

sub _annotate_feature {
    my ($self, $annotation) = @_;
    $annotation->clear_annotation;

    my @fields = PFunc::Annotation::get_valid_fields;

    foreach my $field ( @fields ) {
        my $value;
        if( exists( $transmembrane->{$field} ) ){
            $value = $transmembrane->{$field};
        } else {
            $value = "";
        }

        $annotation->set( $field, $value, $source, $annotation_type );
    }
}


