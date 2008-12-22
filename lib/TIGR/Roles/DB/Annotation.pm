package TIGR::Roles::DB::Annotation;

## Works with schema Chado, vendor Mysql
## could be generalized in the future if need be

use strict;
use warnings;
use DBI;
use Carp;
use Data::Dumper;

sub new {
    my ($class, %args) = @_;
    my $self = {};
    bless($self, $class);
    $self->_init(%args);
    return $self;
}
sub DESTROY {
    my $self = shift;
    my $dbh = $self->_param('dbh');
    $dbh->disconnect();
}

sub _init {
    my ($self, %args) = @_;

    if( $args{'user'} ) {
        $self->user($args{'user'});
    }
    if( $args{'password'} ) {
        $self->password($args{'password'});
    }
    if( $args{'database'} ) {
        $self->database($args{'database'});
    }
    if( $args{'server'} ) {
        $self->server($args{'server'});
    }    

    $self->_connect();
}

sub get_annotation {
    my ($self, $annotation_feature_type) = @_;
    $annotation_feature_type = "transcript" unless( $annotation_feature_type );

    my $query = "SELECT f.uniquename, fp1.value, fp2.value ".
        "FROM feature f, featureprop fp1, featureprop fp2, cvterm cf, cvterm cfp1, cvterm cfp2 ".
        "WHERE f.feature_id = fp1.value ".
        "AND f.type_id = cf.cvterm_id ".
        "AND cf.name = '$annotation_feature_type' ".
        "AND fp1.type_id = cfp1.cvterm_id ".
        "AND cfp1.name = 'gene_product_name' ".
        "AND f.feature_id = fp2.feature_id ".
        "AND fp2.type_id = cfp2.cvterm_id ".
        "AND cfp2.name = 'gene_product_name_source'";

    my $results = $self->_do_select_query( $query );
    my %hash = map { $_->[0] => { 'gene_product_name' => $_->[1],
                                  'gene_product_name_source' => $_->[2] } } @{$results};
    return \%hash;
    
}

sub get_transcript_to_CDS_mapping {
    my ($self) = @_;

    my $query = "SELECT f1.uniquename, c1.name, f2.uniquename ".
        "FROM feature_relationship fr, feature f1, feature f2, cvterm c1, cvterm c2, cvterm c3 ".
        "WHERE fr.subject_id = f1.feature_id ".
        "AND fr.type_id = c3.cvterm_id ".
        "AND f2.feature_id = fr.object_id ".
        "AND c3.name = 'derives_from' ".
        "AND f1.type_id = c1.cvterm_id ".
        "AND f2.type_id = c2.cvterm_id ".
        "AND c1.name = 'CDS' ".
        "AND c2.name = 'transcript'";

    my $results = $self->_do_select_query( $query );
    my %retval = map { $_->[2] => $_->[0] } @{$results};
    return \%retval;
    
}

sub get_transcript_to_polypeptide_mapping {
    my ($self) = @_;

    my $query = "SELECT f1.uniquename, f2.uniquename ".
        "FROM feature f1, feature f2, cvterm c1, cvterm c2, feature_relationship fr, cvterm c3 ".
        "WHERE f1.feature_id = fr.subject_id ".
        "AND f2.feature_id = fr.object_id ".
        "AND fr.type_id = c3.cvterm_id ".
        "AND c3.name = 'part_of' ".
        "AND f1.type_id = c1.cvterm_id ".
        "AND c1.name = 'polypeptide' ".
        "AND f2.type_id = c2.cvterm_id ".
        "AND c2.name = 'transcript'";

    my $results = $self->_do_select_query( $query );
    my %retval = map { $_->[1] => $_->[0] } @{$results};
    return \%retval;
    
}

sub get_all_gene_product_name_source {
    my ($self) = @_;

    my $query = "SELECT f.uniquename, fp.value ".
        "FROM feature f, featureprop fp, cvterm c ".
        "WHERE fp.feature_id = f.feature_id ".
        "AND fp.type_id = c.cvterm_id ".
        "AND c.name = 'gene_product_name_source'";

    my $results = $self->_do_select_query($query);
    my %retval = map { $_->[0] => $_->[1] } @{$results};
    return \%retval;
}

sub get_all_common_names {
    my ($self, $type) = @_;
    
    my $type_cvterm_query = "SELECT cvterm_id ".
        "FROM cvterm ".
        "WHERE name = '$type'";

    my $tq_results = $self->_do_select_query( $type_cvterm_query );
    my $type_cvterm_id = $tq_results->[0]->[0];
    croak("Could not select cvterm_id for name = '$type'") unless( $type_cvterm_id );

    my $query = "SELECT f.uniquename, fp.value ".
        "FROM feature f, featureprop fp, cvterm c ".
        "WHERE c.cvterm_id = fp.type_id ".
        "AND c.name = 'gene_product_name' ".
        "AND fp.feature_id = f.feature_id ".
        "AND f.type_id = $type_cvterm_id";

    my $results = $self->_do_select_query($query);
    my %retval = map { $_->[0] => $_->[1] } @{$results};
    return \%retval;
}

sub _do_insert_query { 
    my ($self, $query) = @_;

    my $dbh = $self->dbh();
    my $sth = $dbh->prepare( $query );

    eval {
        $sth->execute();
        $dbh->commit();
    };
    if( $@ ) {
        $dbh->rollback();
        croak("Database error: $DBI::errstr on query $query\n");
    }

    $sth->finish();
    
}

sub _do_delete_query {
    my ($self, $query) = @_;
    my $num_rows_affected = 0;

    my $dbh = $self->dbh();

    eval {
        $num_rows_affected = $dbh->do($query);
        $dbh->commit();
    };
    if( $@ ) {
        $dbh->rollback();
        croak("Database error: $DBI::errstr on query $query\n");
    }

    return $num_rows_affected;
}

sub _do_select_query {
    my ($self, $query, $column_name) = @_;
    my $dbh = $self->dbh();
    my $sth;
    
    eval {
        $sth = $dbh->prepare( $query );
        $sth->execute();
    };
    if( $@ ) {
        croak("Problem selecting with query [$query] ".
            DBI->errstr);
    }

    my $retval = $sth->fetchall_arrayref();

    $sth->finish;
    
    return $retval;
}

sub remove_tigr_roles {
    my ($self) = @_;
    my $dbh = $self->{'_dbh'};

    #Query out all tigr role annotations
    my $cv_name = "TIGR_role";

    #from the feature_cvterm_prop table
    my $fcp_query = "SELECT fcp.feature_cvtermprop_id ".
        "FROM feature_cvtermprop fcp, feature_cvterm fc, cvterm ct, cv ".
        "WHERE cv.name = 'TIGR_role' ".
        "AND cv.cv_id = ct.cv_id ".
        "AND ct.cvterm_id = fc.cvterm_id ".
        "AND fcp.feature_cvterm_id = fc.feature_cvterm_id";

    my $fcp_results = $self->_do_select_query( $fcp_query );
    my @fcp_ids = map { $_->[0] } @{$fcp_results};

    my $fcp_rows_deleted = 0;
    if( @fcp_ids > 0 ) {
        my $in_clause = join( ', ', @fcp_ids );

        my $fcp_delete = "DELETE FROM feature_cvtermprop ".
            "WHERE feature_cvtermprop_id IN ($in_clause)";
        $fcp_rows_deleted = $self->_do_delete_query( $fcp_delete );
    }
        
    
    #from the feature_cvterm table
    my $fc_query = "SELECT fc.feature_cvterm_id ".
        "FROM feature_cvterm fc, cvterm ct, cv ".
        "WHERE cv.name = 'TIGR_role' ".
        "AND cv.cv_id = ct.cv_id ".
        "AND ct.cvterm_id = fc.cvterm_id";

    my $fc_results = $self->_do_select_query( $fc_query );
    my @fc_ids = map { $_->[0] } @{$fc_results};

    my $fc_rows_deleted = 0;
    if( @fc_ids > 0 ) {
        my $fc_in_clause = join( ', ', @fc_ids );

        my $fc_delete = "DELETE FROM feature_cvterm ".
            "WHERE feature_cvterm_id IN ($fc_in_clause)";
        $fc_rows_deleted = $self->_do_delete_query( $fc_delete );
    }

    return ( $fc_rows_deleted, $fcp_rows_deleted );

}

## Will load an entry in the feature_cvterm table
## as well as the feature_cvtermprop table.
sub annotate_feature_with_role {
    my ($self, $feature_uniquename, $role_id, $ev_accession) = @_;

    croak( "Can not annotate feature with TIGR Role because a parameter is missing ".
           "to the subroutine call. (&annotate_feature_with_role(\$feature_uniqename, ".
           "\$role_id[, \$evidence_accesssion]) )") 
        unless( $feature_uniquename && $role_id );

    my $cvterm_id = $self->get_cvterm_id_by_role_id( $role_id );
    my $feature_id = $self->get_feature_id_by_uniquename( $feature_uniquename );

    return 0 if( $self->does_assignment_exist( $cvterm_id, $feature_id ) );

    #Get the next feature_cvterm_id
    my $feature_cvterm_id_query = "SELECT MAX(feature_cvterm_id) FROM feature_cvterm";
    my $results = $self->_do_select_query( $feature_cvterm_id_query );
    my $next_feature_cvterm_id = $results->[0]->[0] + 1;

    my $query = "INSERT INTO feature_cvterm ( feature_cvterm_id, feature_id, cvterm_id, pub_id ) ".
        "VALUES ( $next_feature_cvterm_id, $feature_id, $cvterm_id, 1)";

    $self->_do_insert_query( $query );


    #If we don't have an ev_accession, our work here is done.
    if( $ev_accession ) {
        #Get the next feature_cvtermprop_id
        my $feature_cvtermprop_id_query = "SELECT MAX(feature_cvtermprop_id) from feature_cvtermprop";
        $results = $self->_do_select_query( $feature_cvtermprop_id_query );
        my $next_feature_cvtermprop_id = $results->[0]->[0] + 1;

        croak("No next_feature_cvtermprop_id") unless( $next_feature_cvtermprop_id );
        croak("No next_feature_cvterm_id") unless( $next_feature_cvterm_id );
        croak("No ev_accession") unless( $ev_accession );
        
        my $fc_query = "INSERT INTO feature_cvtermprop ".
            "( feature_cvtermprop_id, feature_cvterm_id, type_id, value ) ".
            "SELECT $next_feature_cvtermprop_id, $next_feature_cvterm_id, c.cvterm_id, '$ev_accession' ".
            "FROM cvterm c ".
            "WHERE c.name = 'inferred from electronic annotation'";
        
        $self->_do_insert_query( $fc_query );
    }

    return 1;
    
}

sub insert_general_hypo_role {
    my ($self) = @_;

    #get the db_id for tigr roles
    my $tigr_roles_db_id_query = "SELECT db_id from db where name = 'TIGR_role'";
    my $results = $self->_do_select_query( $tigr_roles_db_id_query );
    my $tigr_roles_db_id = $results->[0]->[0];

    #get the cv_id for tigr roles
    my $tigr_roles_cv_id_query = "SELECT cv_id from cv where name = 'TIGR_role'";
    $results = $self->_do_select_query( $tigr_roles_cv_id_query );
    my $tigr_roles_cv_id = $results->[0]->[0] ;
    
    #get the next dbxref_id
    my $dbxref_id_query = "SELECT MAX(dbxref_id) from dbxref";
    $results = $self->_do_select_query( $dbxref_id_query );
    my $next_dbxref_id = $results->[0]->[0] + 1;

    #insert into dbxref if there's not already this entry
    my $dbxref_exist_query = "SELECT dbxref_id from dbxref ".
        "WHERE db_id = $tigr_roles_db_id ".
        "AND accession = '856'";
    $results = $self->_do_select_query( $dbxref_exist_query );
    
    my $insert;
    if( $results->[0]->[0] ) {
        $next_dbxref_id = $results->[0]->[0];
    } else {
        $insert = "INSERT INTO dbxref ( dbxref_id, db_id, accession ) ".
            "VALUES ( $next_dbxref_id, $tigr_roles_db_id, '856' )";
        $self->_do_insert_query( $insert );
    }

    #get the next cvterm_id
    my $cvterm_id_query = "SELECT MAX( cvterm_id ) from cvterm";
    $results = $self->_do_select_query( $cvterm_id_query );
    my $next_cvterm_id = $results->[0]->[0] + 1;

    #insert into cvterm
    my $cvterm_exists_query = "SELECT cvterm_id FROM cvterm ".
        "WHERE cv_id = $tigr_roles_cv_id ".
        "AND name = 'General Hypothetical Protein' ".
        "AND dbxref_id = $next_dbxref_id";
    $results = $self->_do_select_query( $cvterm_exists_query );

    if( $results->[0]->[0] ) {
        $next_cvterm_id = $results->[0]->[0];
    } else {
        $insert = "INSERT INTO cvterm ( cvterm_id, cv_id, name, dbxref_id ) ".
            "VALUES ( $next_cvterm_id, $tigr_roles_cv_id, 'General Hypothetical Protein', $next_dbxref_id )";
        $self->_do_insert_query( $insert );
    } 

    #get the next cvterm_dbxref_id
    my $cvterm_dbxref_id_query = "SELECT MAX( cvterm_dbxref_id ) from cvterm_dbxref";
    $results = $self->_do_select_query( $cvterm_dbxref_id_query );
    my $next_cvterm_dbxref_id = $results->[0]->[0] + 1;

    #insert into cvterm_dbxref if it doesn't already exist
    my $cvterm_dbxref_exists_query = "SELECT cvterm_dbxref_id FROM cvterm_dbxref ".
        "WHERE cvterm_id = $next_cvterm_id ".
        "AND dbxref_id = $next_dbxref_id";
    $results = $self->_do_select_query( $cvterm_dbxref_exists_query );
    
    if( $results->[0]->[0] ) {
        $next_cvterm_dbxref_id = $results->[0]->[0];
    } else {
        $insert = "INSERT INTO cvterm_dbxref ( cvterm_dbxref_id, cvterm_id, dbxref_id ) ".
            "VALUES ( $next_cvterm_dbxref_id, $next_cvterm_id, $next_dbxref_id )";
        $self->_do_insert_query( $insert );
    } 

}

sub does_assignment_exist {
    my ($self, $cvterm_id, $feature_id) = @_;

    croak("subroutine &does_assignemnt_exist requires a cvterm_id passed in")
        unless( $cvterm_id );
    
    my $query = "SELECT feature_cvterm_id ".
        "FROM feature_cvterm ".
        "WHERE feature_id = $feature_id ".
        "AND cvterm_id = $cvterm_id";

    my $results = $self->_do_select_query( $query );
    return scalar( @{$results} );

}

sub get_cvterm_id_by_role_id {
    my ($self, $role_id) = @_;

    my $query = "SELECT c.cvterm_id ".
        "FROM dbxref d1, dbxref d2, db, cvterm_dbxref cd, cvterm c ".
        "WHERE d1.accession = '$role_id' ".
        "AND d1.db_id = db.db_id ".
        "AND db.name = 'TIGR_role' ".
        "AND cd.dbxref_id = d1.dbxref_id ".
        "AND cd.cvterm_id = c.cvterm_id ".
        "AND c.dbxref_id = d2.dbxref_id ".
        "AND d2.db_id = db.db_id";
    
    my $results = $self->_do_select_query( $query );
    my $cvterm_id;
    $cvterm_id = $results->[0]->[0];

    return $cvterm_id;
    
}

sub get_feature_id_by_uniquename {
	my ($self, $uniquename) = @_;
	
	my $query = "SELECT feature_id ".
        "FROM feature ".
        "WHERE uniquename = '$uniquename'";

    my $results = $self->_do_select_query( $query );
    my $feature_id;
    $feature_id = $results->[0]->[0];
    return $feature_id;
	
}

sub dbh {
    my ($self) = @_;
    return $self->_param('dbh');
}

sub user {
    my ($self,$user) = @_;
    return $self->_param('user',$user);
}
sub password {
    my ($self,$password) = @_;
    $self->_param('password',$password);
    return undef;
}
sub database {
    my ($self, $database) = @_;
    return $self->_param('database',$database);
}
sub server {
    my ($self, $server) = @_;
    return $self->_param('server', $server);
}
sub _connect {
    my ($self) = @_;
    my $db = $self->database();
    my $user = $self->user();
    my $password = $self->_param('password');
    my $server = $self->server();

    my $dbh;
    eval {
        $dbh = DBI->connect("DBI:mysql:$db:$server", "$user", "$password",
                            { 
                                'RaiseError' => 1,
                                'AutoCommit' => 0,
                          } );
    };
    if( $@ ) {
        carp("Could not connect to database ".DBI->errstr);
    }

    $self->_param('dbh', $dbh);
}
sub _param {
    my ($self, $name, $value) = @_;
    if( $value ) {
        $self->{"_$name"} = $value;
    } else {
        return $self->{"_$name"};
    }
}

1;
