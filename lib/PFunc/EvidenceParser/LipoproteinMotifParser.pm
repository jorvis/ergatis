package PFunc::EvidenceParser::LipoproteinMotifParser;

use strict;
use warnings;
use XML::Twig;
use Data::Dumper;

use base qw( PFunc::EvidenceParser );

########## Class Vars ##########
my $annotation_type = "LipoproteinMotif";
my $lipoprotein = {
    'gene_product_name' => "putative lipoprotein",
    'GO'                => ['GO:0008150', 'GO:0003674', 'GO:0016020'],
    'TIGR_Role'         => 88,
};
my $source = "LipoproteinMotif";
################################

sub new {
    my ($class, %args) = @_;
    my $self = $class->SUPER::new( %args );
    bless($self, $class);
    #$self->_init_lipop_parser( %args );
    return $self;
}

sub _parse {
    my ($self, $fh) = @_;
    
    my $twig = new XML::Twig( 'twig_roots' => {
        'Sequence[@class="polypeptide"]' => sub { $self->_handle_feature( @_ ) },
    });

    $twig->parse($fh);
}

sub _handle_feature {
    my ($self, $twig, $el) = @_;
    
    my $annote_id = $self->lookup_feature_id( $el->att('id'), $self->annotate_on );
    my $annotation = $self->get_feature_annotation( $annote_id );
    return unless( $annotation );

    my @fields = PFunc::Annotation::get_valid_fields;

    foreach my $field ( @fields ) {
        my $value;
        if( exists( $lipoprotein->{$field} ) ) {
            $value = $lipoprotein->{$field};
        } else {
            $value = "";
        }
        $annotation->set( $field, $value, $source, $annotation_type );
    }

}
