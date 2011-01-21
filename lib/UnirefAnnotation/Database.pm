package UnirefAnnotation::Database;

use strict;
use warnings;
use Data::Dumper;
use UnirefAnnotation::DBUtil;

sub new {
    my ($class, %opts) = @_;
    my $self = {};
    bless($self, $class);
    $self->_init( \%opts );
    return $self;
}

sub get_annotation {
    my ($self, $cluster_acc) = @_;
    my $udb = $self->{'_db'};

    my $clust_id = $udb->get_cluster_id_by_acc( $cluster_acc );
    return [] unless( defined( $clust_id ) );
    
    my @assertions = $udb->get_cluster_assertions( $clust_id );
    return [] unless( keys %assertions );

    return \@assertions;
}

sub is_trusted {
    my ($self, $cluster_acc) = @_;
    my $udb = $self->{'_db'};

    my $clust_id = $udb->get_cluster_id_by_acc( $cluster_acc );
    return 0 unless( defined( $clust_id ) );
    
    return $udb->cluster_is_trusted( $clust_id );
}

sub _init {
    my ($self, $opts) = @_;

    #required:: username password
    #optioonal: test
    die("Option username is  required") unless( exists( $opts->{'username'} ) );
    die("Option password is  required") unless( exists( $opts->{'password'} ) );
    my $test = 0;
    $test = $opts->{'test'} if( exists( $opts->{'test'} ) );

    my $db = new UnirefAnnotation::DBUtil( 'username' => $opts->{'username'},
                                           'password' => $opts->{'password'},
                                           'testt' => $test );
    $self->{'_db'} = $db;
                                                 
}

1;
