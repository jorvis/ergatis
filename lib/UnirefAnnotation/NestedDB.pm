package UnirefAnnotation::NestedDB;

use strict;
use warnings;
use NestedDB;
use Data::Dumper;

sub new {
    my ($class, %opts) = @_;
    my $self = {};
    bless($self, $class);
    $self->_init( \%opts );
    return $self;
}

sub get_annotation {
    my ($self, $cluster_acc) = @_;
    my $db = $self->{'_db'};

    return [] unless( exists( $db->{$cluster_acc} ) );
    return $db->{$cluster_acc}->{'assertions'};
}

sub is_trusted {
    my ($self, $cluster_acc) = @_;
    my $db = $self->{'_db'};
    return 0 unless( exists( $db->{$cluster_acc} ) );
    return $db->{$cluster_acc}->{'is_trusted'};
}

sub _init {
    my ($self, $opts) = @_;

    #required:: path
    die("Path to the lookup file is required")
        unless( exists( $opts->{'path'} ));
    
    my %db;
    tie( %db, "NestedDB", $opts->{'path'} ) or die("Couldn't tie hash to file $!");
    
    $self->{'_db'} = \%db;
}

1;
