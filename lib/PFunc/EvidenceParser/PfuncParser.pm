package PFunc::EvidenceParser::PfuncParser;

use strict;
use warnings;
use Carp;
use XML::Twig;
use File::OpenFile qw( open_file );
use PFunc::Annotation;
use Data::Dumper;

use base qw(PFunc::EvidenceParser);

### class variables ###
my $annotation_type = "Pfunc";

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
        "Feature" => sub { $self->_feature_handler(@_) },
    });

    $twig->parse( $fh );

}

sub _feature_handler {
    my ($self, $twig, $feature_elem) = @_;

    my $feat_id = $feature_elem->att('id');
    my $class = $feature_elem->att('class');

    ## make sure we annotate on the correct feature type
    my $annotate_id;
    eval {
        $annotate_id = $self->lookup_feature_id( $feat_id, $self->annotate_on );
    };
    if( $@ ) {
        croak("$@");
    }
    my $annotation = $self->get_feature_annotation( $annotate_id );
    
    # get all children Attribute elements
    my @atts = $feature_elem->children( 'Attribute' );
    
    foreach my $attribute ( @atts ) {
        $self->_handle_attribute( $annotation, $attribute );
    }
    
    # get all children Attribute-list elements
    my @att_lists = $feature_elem->children( 'Attribute-list' );

    foreach my $attribute_list ( @att_lists ) {
        $self->_handle_attribute_list( $annotation, $attribute_list );
    }

}

sub _handle_attribute {
    my ($self, $annotation, $attribute) = @_;

    my $name = $attribute->att('name');
    my $content = $attribute->att('content');

    if( $name =~ /(.*)_source/ ) {
        my $field = $1;
        print "Calling set on $field for source\n";
        $annotation->set( $field, undef, $content, $annotation_type );
    } else {
        print "Calling set on $name\n";
        $annotation->set( $name, $content, undef, $annotation_type, undef );
    }
}

sub _handle_attribute_list {
    my ($self, $annotation, $attribute_list) = @_;
    
    my @atts = $attribute_list->children( 'Attribute' );

    # the first attribute should define the annotation type
    my $first_att = $atts[0];
    my $annotation_name = $first_att->att('name');
    my $annotation_content = $first_att->att('content');
    
    # the second atribute should contain the source and ev_code
    my $sec_att = $atts[1];
    my $ev_code = $sec_att->att('name');
    my $source = $sec_att->att('content');

    # add the annotation
    $annotation->add( $annotation_name, $annotation_content, $source, $annotation_type );

}

1;
