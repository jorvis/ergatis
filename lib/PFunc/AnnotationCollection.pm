package PFunc::AnnotationCollection;

use strict;
use warnings;
use Data::Dumper;
use Carp;
use File::OpenFile qw( open_file );
use PFunc::Annotation;
use PFunc::AnnotationOrder qw(@annotation_order);

our @annotation_order;
my $default_annotation_order = "AnnotationOrder.pm";
my %max_ranks;

sub new {
    my ($class, %params) = @_;
    my $self = {};
    bless($self, $class);
    $self->_init( %params );
    return $self;
}

sub DESTROY {
    my ($self) = @_;
    if( $self->{'_fh'} ) {
        close( $self->{'_fh'} );
    }
}

sub annotations_from_file {
    my ($self, $file) = @_;
    $self->{'_fh'} = open_file( $file, 'in' );
}

sub next_annotation_from_file {
    my ($self) = @_;
    my $fh = $self->{'_fh'};
    my ($retval, $string);
    $retval = PFunc::Annotation::to_object( $string ) 
        if( defined( $string = <$fh> ) );
    return $retval;
}

sub annotation_order {
    my ($self, $ds) = @_;
    if( defined($ds) ) {
        $self->{'_annotation_order'} = $ds;
    }
    return $self->{'_annotation_order'};
}

sub add_annotation {
    my ($self, $annotation) = @_;
    my $final_annot = $annotation;
    if( $self->get_annotation( $annotation->get_feature_id ) ) {
        $final_annot = $self->_merge_annotations( $annotation, 
                                                  $self->get_annotation( $annotation->get_feature_id ) );
    }
    $self->{"_annotations"}->{$annotation->get_feature_id} = $final_annot;
    return $final_annot;
}

sub get_annotation {
    my ($self, $feature_id) = @_;
    return ( exists($self->{'_annotations'}->{$feature_id}) ) ? 
        $self->{'_annotations'}->{$feature_id} :
        undef;
}

sub annotations {
    my ($self) = @_;
    return values( %{$self->{'_annotations'}} );
}

sub output {
    my ($self, $ofh) = @_;
    if( defined( $ofh ) ) {
        $self->{'_output'} = $ofh;
    }
    return $self->{'_output'};
}

sub flush {
    my ($self) = @_;
    my $ofh = $self->output;

    foreach my $annot ( $self->annotations ) {
        print $ofh $annot->to_string()."\n";
        $self->{'_total_annotations_printed'}++;
    }

    $self->{'_annotations'} = {};
}

sub total_annotations_printed {
    my ($self) = @_;
    return $self->{'_total_annotations_printed'};
}

sub annotations_in_buffer {
    my ($self) = @_;
    return scalar( keys %{$self->{'_annotations'}} );
}

sub _merge_annotations {
    my ($self, $new_annotation, $current_annotation, $feature_id) = @_;

    #If they didn't pass in a feature id, then use the feature id from the first
    #annotation.  If that doesn't exist (which it always should) it will take
    #from the second annotation.  These should usually be the same.
    if( !defined( $feature_id ) ) {
        $feature_id = ( $current_annotation->get_feature_id ) ? 
            $current_annotation->get_feature_id : $new_annotation->get_feature_id;
    }

    my $ret_annot = new PFunc::Annotation( 'feature_id' => $feature_id );

    my @fields = PFunc::Annotation::get_sorted_fields();
    
    foreach my $field ( @fields ) {
        
        #If they both have the annotation
        if( $new_annotation->is_annotated( $field ) && $current_annotation->is_annotated( $field ) ) {
            my $cur_rank = $self->_get_annot_rank( $field, $current_annotation->_get_type( $field ) );
            my $new_rank = $self->_get_annot_rank( $field, $new_annotation->_get_type( $field ) );

            if( $cur_rank <= $new_rank ) {
                $ret_annot->set( $field, $current_annotation->get( $field ) );
            } elsif( $new_rank < $cur_rank ) {
                $ret_annot->set( $field, $new_annotation->get( $field ) );
            } else {
                confess("Could not merge annotations");
            }
        } elsif( $new_annotation->is_annotated($field) ) {
            $ret_annot->set( $field, $new_annotation->get( $field ) );
        } elsif( $current_annotation->is_annotated($field) ) {
            $ret_annot->set( $field, $current_annotation->get($field) );
        }

    }

    unless( $ret_annot->has_annotation ) {
        print Dumper( $ret_annot );
        croak("At the end of merge, retval has no annotations");
    }

    return $ret_annot;
}

sub _get_annot_full_or_partial {
    my ($self, $type) = @_;
    my $retval;

    my $ao = $self->annotation_order();
    
    if( exists( $ao->{$type} ) ) {
        $retval = $ao->{$type}->{'_full_or_partial'};
    } else {
        confess("Annotation type $type cannot be found in the AnnotationOrder");
    }

    return $retval;
}

sub _get_annot_rank {
    my ($self, $field, $type) = @_;
    my $retval;

    my $ao = $self->annotation_order();

    if( exists( $ao->{$type} ) ) {
        confess("The field $field was not defined in AnnotationOrder for annotation type $type")
            unless( $retval = $ao->{$type}->{'_fields'}->{$field} );
    } else {
        $retval = $max_ranks{ $field };
    }
    
    return $retval;
    
}

#Will give weights to the ordered annotation types
sub _reformat_annotation_order {
    my ($ordered) = @_;
    my $retval;
    
    my @fields = PFunc::Annotation::get_valid_fields();
    my %indices;
    map { $indices{$_} = 1 } @fields;

    foreach my $annotRank ( @{$ordered} ) {
        my ($type, $ranked_fields, $full_or_partial) = @{$annotRank};

        foreach my $ranked_field( @{$ranked_fields} ) {

            if( $ranked_field eq 'all' ) {
                map { $retval->{ $type }->{'_fields'}->{$_} = $indices{$_}++
                          unless( exists( $retval->{$type}->{'_fields'}->{$_} ) );
                      } @fields;
                $retval->{$type}->{'_all'} = 1;
                
                last;
            }
            confess("Could not find field $ranked_field in Annotation definition.  This was found in ".
                    "AnnotationOrder for type $type") unless( $indices{$ranked_field} );
            $retval->{$type}->{'_fields'}->{$ranked_field} = $indices{$ranked_field}++;

        }
    }

    foreach my $key ( keys %indices ) {
        $max_ranks{ $key } = $indices{ $key };
    }
    
    return $retval;
    
}

sub _init {
    my ($self, %params) = @_;
    
    if( $params{'annotation_order'} ) {
        $self->annotation_order( &_reformat_annotation_order( $params{'annotation_order'} ) );
    } else {
        $self->annotation_order( &_reformat_annotation_order( \@annotation_order ) );
    }

    if( $params{'output'} ) {
        $self->output( $params{'output'} );
    } else {
        $self->output( *STDOUT );
    }

    $self->{'_annotations'} = {};
    $self->{'_total_annotations_printed'} = 0;
}

1;
