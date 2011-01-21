package UnirefAnnotation::DBUtil;

use strict;
use warnings;
use DBI;

my $DB = "uniprot_annotation";
my $TEST_DB = "uniprot_annotation_test";
my $HOST = "jabba";

sub new {
    my ($class, %opts) = @_;
    my $self = {};
    bless( $self, $class );
    $self->_init( \%opts );
    return $self;
}

############################################################
#                Generic Subs                              #
############################################################
sub do_select_query {
    my ($self, $query, @args) = @_;
    my $dbh = $self->_param('dbh');
    my $sth = $dbh->prepare( $query );

    $sth->execute(@args) or
        die("Could not execute query [$query]: ".DBI->errstr);

    my $results = $sth->fetchall_arrayref;
    $sth->finish();
    return $results;
}

sub do_single_select {
    my ($self, $query, @args) = @_;
    my $results = $self->do_select_query( $query, @args );
    my $retval;
    if( @{$results} != 0 ) {
        $retval = $results->[0]->[0];
    }
    return $retval;
}

sub do_insert_query {
    my ($self, $query, @args) = @_;
    my $dbh = $self->_param('dbh');
    my $sth = $dbh->prepare( $query );
    
    $sth->execute( @args ) or
        die("Could not execute query [$query]: ".DBI->errstr);
    $sth->finish();

    $dbh->{'mysql_insertid'}; 
}

sub do_delete_or_update_query {
    my ($self, $query, @args) = @_;
    my $dbh = $self->_param('dbh');
    my $sth = $dbh->prepare( $query );
    
    $sth->execute( @args ) or
        die("Could not execute query [$query]: ".DBI->errstr);
    $sth->finish();
}

sub do_update_query {
    my ($self, $query, @args) = @_;
    $self->do_delete_or_update_query( $query, @args );
}

sub do_delete_query {
    my ($self, $query, @args) = @_;
    $self->do_delete_or_update_query( $query, @args );
}
############################################################
#                Select Methods                            #
############################################################
sub get_cluster_type_id {
    my ($self, $cluster_type) = @_;
    my $query = "SELECT id FROM uniref_cluster_types ".
        "WHERE name = ?";
    $self->do_single_select( $query, $cluster_type );
}

sub does_cluster_type_id_exist{
    my ($self, $cluster_type_id) = @_;
    my $query = "SELECT id FROM uniref_cluster_types WHERE id = ?";
    $self->do_single_select( $query, $cluster_type_id );
}

sub does_cluster_id_exist {
    my ($self, $cluster_id) = @_;
    my $query = "SELECT id FROM uniref_clusters WHERE id = ?";
    my $retval = $self->do_single_select( $query, $cluster_id );
    $retval;
}

sub get_cluster_member_id_by_accession {
    my ($self, $acc) = @_;
    my $query = "SELECT id FROM uniref_cluster_members ".
        "WHERE accession = ?";
    $self->do_single_select( $query, $acc );
}

sub get_cluster_id_by_acc {
    my ($self, $cluster_acc) = @_;
    my $query = 
        "SELECT id FROM uniref_clusters ".
        "WHERE accession = ?";
    $self->do_single_select( $query, $cluster_acc );
}

sub get_cluster_id_by_member_acc {
    my ($self, $cluster_member_acc) = @_;
    my $query = 
        "SELECT cluster_id FROM uniref_cluster_members ".
        "WHERE accession = ?";
    $self->do_single_select( $query, $cluster_member_acc );
}

sub get_cluster_assertions {
    my ($self, $cluster_id) = @_;
    my $query = 
        "SELECT type, value, is_manual, source, assigned_by, id FROM assertions ".
        "WHERE cluster_id = ?";
    my $results = $self->do_select_query( $query, $cluster_id );

    my @retval;
    foreach my $res ( @{$results} ) {
        my $tmp = {
			'type' => $res->[0],
            'value' => $res->[1],
            'is_manual' => $res->[2],
            'source' => $res->[3],
            'assigned_by' => $res->[4],
            'id' => $res->[5]
            };
		push(@retval, $tmp);
    }
    return @retval;
}

sub get_cluster_assertions_by_type {
    my ($self, $cluster_id, $type) = @_;
    my %assertions = $self->get_cluster_assertions( $cluster_id );
    return unless( exists( $assertions{$type} ) );
    return $assertions{$type};
}

sub cluster_is_trusted {
    my ($self, $cluster_id) = @_;
    my $query = 
        "SELECT is_trusted FROM uniref_clusters ".
        "WHERE id = ?";
    my $is_trusted = $self->do_single_select( $query, $cluster_id );
    $is_trusted;
}

sub get_cluster_members {
    my ($self, $cluster_id) = @_;
    my $query = 
        "SELECT accession, id FROM uniref_cluster_members ".
        "WHERE cluster_id = ?";
    my $results = $self->do_select_query( $query, $cluster_id );
    my @members;
    map { push(@members, $_->[0] ); } @{$results};
    return @members;
}

sub get_all_clusters_with_assertions {
    my ($self) = @_;
    my $query = 
        "SELECT c.id, c.accession, c.is_trusted, a.type, a.value ".
        "FROM uniref_clusters c, assertions a ".
        "WHERE a.cluster_id = c.id";
    $self->do_select_query( $query );    
}
############################################################
#                      Insert Subs                         #
############################################################
sub insert_cluster_type {
    my ($self, $cluster_type) = @_;
    my $query = "INSERT INTO uniref_cluster_types (name) ".
        "VALUES( ? )";
    $self->do_insert_query( $query, $cluster_type );
}

sub insert_cluster {
    my ($self, $accession, $is_trusted, $cluster_type_id) = @_;
    die("Cluster type id: $cluster_type_id does not exist")
        unless( $self->does_cluster_type_id_exist( $cluster_type_id) );
    my $query = 
        "INSERT INTO uniref_clusters( accession, is_trusted, cluster_type ) ".
        "VALUES( ?, ?, ?)";
    $self->do_insert_query( $query, $accession, $is_trusted, $cluster_type_id );
}

sub insert_cluster_member {
    my ($self, $accession, $cluster_id) = @_;
    die("Cluster id: $cluster_id id does not exist in the database")
        unless( $self->does_cluster_id_exist( $cluster_id ) );
    my $query = 
        "INSERT INTO uniref_cluster_members( accession, cluster_id ) ".
        "VALUES( ?, ?)";
    $self->do_insert_query( $query, $accession, $cluster_id);
}

## Add assign_by sometime
sub insert_assertion {
    my ($self, $cluster_id, $type, $value, $is_manual, $source, $assigned_by) = @_;
    die( "Cluster id: $cluster_id does not exist") 
        unless( $self->does_cluster_id_exist($cluster_id) );
    my $query = 
        "INSERT INTO assertions (cluster_id, type, value, is_manual, source, assigned_by) ".
        "VALUES(?, ?, ?, ?, ?, ?)";
    $self->do_insert_query($query, $cluster_id, $type, $value, $is_manual, $source, $assigned_by);
}

sub set_is_trusted {
    my ($self, $cluster_id, $is_trusted) = @_;
    my $query = "UPDATE uniref_clusters SET is_trusted = ? ".
        "WHERE id = ?";
    $self->do_update_query( $query, $is_trusted, $cluster_id );
}

############################################################
#                    Delete Subs                           #
############################################################
sub remove_assertions_from_cluster {
    my ($self, $cluster_id) = @_;

    my $query = 
        "DELETE FROM assertions ".
        "WHERE cluster_id = ?";
    $self->do_delete_query( $query, $cluster_id );
}

sub remove_assertions_from_cluster_by_type {
    my ($self, $cluster_id, $type) = @_;
    my $query = 
        "DELETE FROM assertions ".
        "WHERE cluster_id = ? AND type = ?";
    $self->do_delete_query( $query, $cluster_id, $type );
}

sub remove_cluster_member_by_acc {
    my ($self, $acc) = @_;
    my $query = 
        "DELETE FROM uniref_cluster_members ".
        "WHERE accession = ?";
    $self->do_delete_query( $query, $acc );
}

############################################################
#                   "Private" Subs                         #
############################################################
sub _param {
    my ($self, $param, $value) = @_;
    if( defined( $value ) ) {
        $self->{$param} = $value;
    }
    return $self->{$param};
}

sub _connect {
    my ($self) = @_;
    my ($db, $host,$user, $pass) = ($self->_param('database'), $self->_param('host'), 
                                    $self->_param('username'), $self->_param('password') );
    my $retval = DBI->connect("dbi:mysql:host=$host;packetSize=8092", $user, $pass, 
                              { 
                                  'AutoCommit' => 1,
                                  'RaiseError' => 1
                                  }
                                      ) or
        die("Could not connect to server:".DBI->errstr);

    $self->_param('dbh', $retval);
    if( defined( $db ) ) {
        $retval->do( "use $db" );
    }
}

sub _init {
    my ($self, $options) = @_;
    
    foreach my $req( qw(username password) ) {
        die("Option $req is required") unless( $options->{$req} );
    }

    my $db = $options->{'database'} || $DB;
    my $host = $options->{'host'} || $HOST;

    $db = $TEST_DB if( $options->{'test'} );
    
    $self->_param('database', $db );
    $self->_param('host', $host );
    $self->_param('username', $options->{'username'} );
    $self->_param('password', $options->{'password'} );
    
    $self->_connect();

    
}

1;
