package Ergatis::Pipeline;

=head1 NAME

Ergatis::Pipeline.pm - A class for representing pipelines

=head1 SYNOPSIS

    use Ergatis::Pipeline;

    my $conf = Ergatis::Pipeline->new( id => '1234', path => '/path/to/some/pipeline.xml' );

=head1 DESCRIPTION

=head1 METHODS

=over 3

=item I<PACKAGE>->new( path => I<'/path/to/some.xml'> )

Returns a newly created "Ergatis::Pipeline" object.

=item I<$OBJ>->path( )

Returns the path to a given pipeline.

=item I<$OBJ>->id( )

Returns the ID of a given pipeline.

=back

=head1 AUTHOR

    Joshua Orvis
    jorvis@tigr.org

=cut

use strict;
use Carp;

umask(0000);

## class data and methods
{

    my %_attributes = (
                        path  => undef,
                        id    => undef,
                      );

    ## class variables
    ## currently none

    sub new {
        my ($class, %args) = @_;
        
        ## create the object
        my $self = bless { %_attributes }, $class;
        
        ## set any attributes passed, checking to make sure they
        ##  were all valid
        for (keys %args) {
            if (exists $_attributes{$_}) {
                $self->{$_} = $args{$_};
            } else {
                croak("$_ is not a recognized attribute");
            }
        }
        
        return $self;
    }
    
    sub run {
        ## path must be defined and exist to run
        if (! defined $_[0]->{path} ) {
            croak("failed to run pipeline: no path set yet.");
        } elsif (! -f $_[0]->{path} ) {
            croak("failed to run pipeline: pipeline file $_[0]->{path} doesn't exist");
        }
        
        
    }
    
    ## accessors
    sub id { return $_[0]->{id} }
    sub path { return $_[0]->{path} }
    
    ## modifiers
    sub set_id { $_[0]->{id} = $_[1] }
    sub set_path { $_[0]->{path} = $_[1] }
}

1==1;

