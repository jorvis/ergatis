package PFunc::EvidenceParser::PriamECParser;

use strict;
use warnings;
use Carp;
use XML::Twig;
use File::OpenFile qw( open_file );
use PFunc::Annotation;
use Data::Dumper;

use base qw(PFunc::EvidenceParser);

### class variables ###
my $annotation_type = "Priam";

sub new {
    my ($class, %args) = @_;

    my $self = $class->SUPER::new( %args );
    bless($self, $class);
    return $self;
}

##### private subroutines #####
sub _parse {
    my ($self, $fh) = @_;
    
    my $twig = new XML::Twig( 'twig_handlers' => {
        "Sequence" => sub { $self->_sequence_handler(@_) },
    });

    $twig->parse( $fh );

}

#only looking for ec numbers
sub _sequence_handler {
    my ($self, $twig, $sequence_elem) = @_;

    my $seq_id = $sequence_elem->att('id');    

    # add the type to the annotation.  This is used when
    # resolving conflicts.
    my $annotation_id = $self->lookup_feature_id( $seq_id, $self->annotate_on );
    my $annotation = $self->get_feature_annotation( $annotation_id );

    # get all children Attribute-list elements
    my @att_lists = $sequence_elem->children( 'Attribute-list' );

    foreach my $attribute_list ( @att_lists ) {
        $self->_handle_attribute_list( $annotation_id, $attribute_list );
    }

}

sub _handle_attribute_list {
    my ($self, $seq_id, $attribute_list) = @_;

    my $annotation = $self->get_feature_annotation( $seq_id );
    
    my @atts = $attribute_list->children( 'Attribute' );

    # the first attribute should define the annotation type
    my $first_att = $atts[0];
    my $annotation_name = $first_att->att('name');
    return unless( $annotation_name eq 'EC' );       # only ec numbers

    my $annotation_content = $first_att->att('content');
    
    # the second atribute should contain the source and ev_code
    my $sec_att = $atts[1];
    my $ev_code = $sec_att->att('name');
    my $source = $sec_att->att('content');

    # add the content, the source and the source type (i.e. priam ec assignment)
    $annotation->set( $annotation_name, $annotation_content, $source, $annotation_type );

}

1;
