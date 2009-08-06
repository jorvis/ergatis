package PFunc::EvidenceParser::Hypothetical;

use strict;
use warnings;
use XML::Twig;
use base qw( PFunc::EvidenceParser );

##### class vars #####
my $annotation_type = "hypothetical";
my $hypo = {
    'gene_product_name' => "hypothetical protein",
    'GO'                => ['GO:0008150', 'GO:0003674', 'GO:0005575'],
    'TIGR_Role'         => 856,
};
my $source = "default";
##########

sub new {
    my ($class, %args) = @_;
    my $self = $class->SUPER::new( %args );
    bless($self, $class);
    return $self;
}

sub _parse {
    my ($self, $fh) = @_;
    
    my $twig = new XML::Twig( 'twig_roots' => {
        'Feature[@class="gene"]' => sub {
            my ($t, $feat) = @_;
            my $annot_id = $self->lookup_feature_id( $feat->att('id'), $self->annotate_on );
            my $annotation = $self->get_feature_annotation( $annot_id );
            
            $self->_set_annotation( $annotation );
            
        }
    });

    $twig->parse( $fh );
}

sub _set_annotation {
    my ($annotation) = @_;

    if( $annotation ) {
        $annotation->clear_annotation;

        my @fields = PFunc::Annotation::get_valid_fields;
        foreach my $field ( @fields ) { 
            my $value = "";
            $value = $hypo->{$field} if( exists( $hypo->{$field} ) );
            $annotation->set( $field, $value, $source, $annotation_type );
        }
        
    }
}

1;
