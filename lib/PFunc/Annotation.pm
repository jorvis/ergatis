package PFunc::Annotation;

# generic implementation of an annotation object
# To add fields, you just need to add it to the %valid_object_variables hash

use strict;
use warnings;
use Carp;
use Class::Struct;
use Data::Dumper;
our $AUTOLOAD;


# define variable types
my ($SINGLE, $MULTIPLE) = (1,2);
my $multi_value_separator = ",";
my $field_separator = "\t";
my $value_source_type_separator = '-;-';

# valid object variables
# defines if the value can have a single value 
# or multiple values and also the order in which they
# will be printed.  I wanted this to be hash, thus
# this weird index thing.
my $index = 0;
my %valid_object_variables = 
    (
     'gene_product_name'             => $index++ ,
     'gene_symbol'                   => $index++ ,
     'EC'                            => $index++ ,
     'GO'                            => $index++ ,
     'TIGR_Role'                     => $index++ ,
     );

sub new {
    my ($class, %args) = @_;
    my $self = {};
    bless( $self, $class );
    $self->_init(%args);
    return $self;
}

sub set {
    my ($self, $field, $value, $source, $type ) = @_;
    return $self->_set_field( $field, $value, $source, $type );
}

sub add {
    my ($self, $field, $value, $source, $type) = @_;
    croak("$field is not a valid field") 
        unless( &is_valid_field( $field ) );
    $self->_add_value( $field, $value );
    $self->_set_source( $field, $source ) if( defined( $source ) );
    $self->_set_type( $field, $type ) if( defined( $type ) );
}

sub get {
    my ($self, $field) = @_;
    return $self->_get_field( $field );
}

sub get_feature_id {
    my ($self) = @_;
    my $retval;
    $retval = $self->{"_feature_id"} if( $self->{"_feature_id"} );
    return $retval;
}
sub set_feature_id {
    my ($self, $feature_id) = @_;
    $self->{"_feature_id"} = $feature_id;
    return $self->{"_feature_id"};
}

sub clear_annotation {
    my ($self) = @_;
    my @fields = &get_valid_fields;
    foreach my $field ( @fields ) {
        $self->_clear_field( $field );
    }
}

sub transfer_annotation {
    my ($self, $annotation) = @_;
    
    $self->set_feature_id( $annotation->get_feature_id );

    foreach my $field( &get_valid_fields ) {
        $self->set( $field, $annotation->get( $field ) );
    }

}

sub to_string {
    my ($self) = @_;

    my @sorted_keys = &get_sorted_fields();

    my @ordered_values = ( $self->get_feature_id );
    my @get = $self->get($sorted_keys[0]);
    
    map { push( @ordered_values, &_value_to_string( $_, $self->get( $_ ) ) || "" ) } @sorted_keys;
    return join ($field_separator, @ordered_values );
}

sub to_object {
    my ($string) = @_;

    my @cols = split( $field_separator, $string, -1 );
    croak("Data string does not correctly represent Annotation object. ".
          "Expecting ".($index+1)." cols and string has ".scalar(@cols).". [$string]")
        unless( @cols == $index+1 );

    my $self = new PFunc::Annotation( 'feature_id' => shift @cols );

    my @sorted_keys = &get_sorted_fields();

    for( my $i = 0; $i < $index; $i++ ) {
        if( defined( $cols[$i] ) ) {
            my $field = $sorted_keys[$i];
            $self->_set_field( $field, &_string_to_value_source_type( $cols[$i] ) );
        }
    }

    return $self;
}

# If this object doesn't really have any information
# i.e if it has no info for anything except feature_id and
# annotation type
sub has_annotation {
    my ($self, @fields) = @_;
    my $has_annotation = 0;

    @fields = &get_valid_fields if( @fields == 0 );

    # cycle through all the valid fields
    foreach my $field ( @fields ) {
        my ($value) = $self->get( $field );
        if( @{$value} > 0 ) {
            $has_annotation = 1;
            last;
        }
    }

    # and return
    return $has_annotation;    
}

# If a field currently has anything for either the value, source or type.
# Takes only one field and it is required.
sub is_annotated {
    my ($self, $field) = @_;
    
    confess( "Field $field is not valid" ) 
        unless( &is_valid_field( $field ) );

    my $is_annotated = 0;

    my ($value, $source, $type) = $self->get($field);
    if( &is_value_defined( $value ) || &is_value_defined( $source ) || 
        &is_value_defined( $type ) ) {
        $is_annotated = 1;
    }

    return $is_annotated;
}

############################## class subroutines ##############################
sub get_sorted_fields {
    return &_get_sorted_fields();
}
sub get_valid_fields {
    return &_get_valid_fields();
}
sub is_valid_field {
    my ($field) = @_;
    my $retval = 0;
    $retval = 1 if( exists( $valid_object_variables{ $field } ) );
    return $retval;
}

############## private subroutines ################
sub _get_sorted_fields {
    return sort { $valid_object_variables{$a} <=> $valid_object_variables{$b} } &get_valid_fields;
}

sub _get_valid_fields {
    return keys %valid_object_variables;
}

sub _string_to_value_source_type {
    my ($string) = @_;
    my ($value, $source, $type) = split( $value_source_type_separator, $string );
    my @ret_value = split( $multi_value_separator, $value );
    return (\@ret_value, $source, $type);
}

sub _value_to_string {
    my ($field, $array_ref, $source, $type) = @_;
    my $retval = "";
    if( defined( $array_ref ) ) {
        my $tmp_string = join( $multi_value_separator, @{$array_ref} );
        $retval = join( $value_source_type_separator, ( $tmp_string, $source, $type ) );
    }
    return $retval;
}

sub _init {
    my ($self, %args) = @_;    
    if( $args{'feature_id'} ) {
        $self->set_feature_id( $args{'feature_id'} );
    }
}

sub _set_field {
    my ($self, $field, $value, $source, $type ) = @_;
    croak("$field is not a valid field") 
        unless( &is_valid_field( $field ) );

    $self->_set_value(  $field, $value  ) if( defined( $value  ) );
    $self->_set_source( $field, $source ) if( defined( $source ) );
    $self->_set_type(   $field, $type   ) if( defined( $type   ) );
    return $value;
}

sub _set_value {
    my ($self, $field, $value) = @_;
    my @push;
    if( ref( $value ) eq 'ARRAY' ) {
        push(@push, @{$value});
    } else {
        push(@push, $value);
    }
    $self->{"_$field"}->{'value'} = \@push;
}

sub _set_source {
    my ($self, $field, $source) = @_;
    $self->{"_$field"}->{'source'} = $source;
}

sub _set_type {
    my ($self, $field, $type) = @_;
    $self->{"_$field"}->{'type'} = $type;
}

sub _get_field {
    my ($self, $field) = @_;
    croak( "Field $field is not a valid parameter for an annotation" )
        unless( &is_valid_field( $field ) );
    return ($self->_get_value($field), $self->_get_source($field), $self->_get_type($field));
}

sub _get_value {
    my ($self, $field) = @_;
    my $retval = [];
    if( exists( $self->{"_$field"} ) && $self->{"_$field"}->{'value'} ) {
        $retval = $self->{"_$field"}->{'value'};
    }
    return $retval;
}

sub _get_source {
    my ($self, $field) = @_;
    my $retval = "";
    if( exists( $self->{"_$field"} ) && $self->{"_$field"}->{'source'} ) {
        $retval = $self->{"_$field"}->{'source'};
    }
    return $retval;
}
sub _get_type {
    my ($self, $field) = @_;
    my $retval = "";
    if( exists( $self->{"_$field"} ) && $self->{"_$field"}->{'type'} ) {
        $retval = $self->{"_$field"}->{'type'};
    }
    return $retval;
}

sub _clear_field {
    my ($self, $field) = @_;
    $self->_clear_value( $field );
    $self->_clear_source( $field );
    $self->_clear_type( $field );
}

sub _clear_value {
    my ($self, $field) = @_;
    $self->{"_$field"}->{'value'} = undef;
}
sub _clear_source {
    my ($self, $field) = @_;
    $self->{"_$field"}->{'source'} = undef;
}
sub _clear_type {
    my ($self, $field) = @_;
    $self->{"_$field"}->{'type'} = undef;
}

sub _add_value {
    my ($self, $field, $value) = @_;
    my @push = ();
    if( ref( $value ) eq 'ARRAY' ) {
        push(@push, @{$value});
    } else {
        push(@push, $value);
    }
    push( @{$self->{"_$field"}->{'value'}}, @push );
    return @push
}
sub is_value_defined {
    my ($value) = @_;
    my $retval = 0;
    if( defined($value) ) {
        if( ref($value) eq 'ARRAY' ) {
            $retval = 1 if( @{$value} > 0 );
        } elsif( ref($value) eq 'HASH' ) {
            $retval = 1 if( keys %{$value} > 0 );
        } elsif( $value && $value !~ /^\s*$/ ) {
            $retval = 1;
        }
    }

    return $retval;
}

sub AUTOLOAD {
    my ($self, @params) = @_;
    
    # grab the name of the call
    my $name = $AUTOLOAD;
    return if( $name =~ /DESTROY/ );
    $name =~ s/.*://;

    # accessor or setter?
    my ($access, $field) = ("","");
    ($access, $field) = ($1,$2) if( $name =~ /^(set|get|add)_(.*)/ );

    # check to make sure it's valid
    my $type;
    confess("Cannot locate object method $name via package \"".__PACKAGE__."\"")
        unless( $type = &is_valid_field( $field ) );

    # determine which method to call
    eval("\$self->_${access}_field( \$field, \@params )");
}
1;
