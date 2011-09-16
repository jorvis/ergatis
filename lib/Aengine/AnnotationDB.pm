package Aengine::AnnotationDB;

use strict;
use warnings;
use DBI;

my @rna_types = qw(rRNA tRNA snRNA ncRNA);
my @molecule_types = qw(assembly pseudomolecule contig chromosome plasmid);

sub new {
    my ($class, %args) = @_;
    my $self = {};
    bless( $self, $class );
    $self->_init(%args);
    return $self;    
}

sub DESTROY {
    my ($self) = @_;
    my $dbh = $self->_param('dbh');
    $dbh->disconnect() if( defined( $dbh ) );
}

###########################################################################
#                           Generic Subroutines                           #
###########################################################################

sub do_select_query {
    my ($self, $query, @params) = @_;
    my $dbh = $self->_param('dbh');
    my $sth = $dbh->prepare( $query );

    eval {
        $sth->execute(@params);
    };
    if( $@ ) {
        die("Could not execute query [$query]: ".DBI->errstr);
    }
    my $results = $sth->fetchall_arrayref;
    $sth->finish();
    return $results;
}

sub change_db {
    my ($self, $db) = @_;
    $self->_param("database", $db);
    my $dbh = $self->_param('dbh');
    $dbh->do("use $db") or
        return 0;
    return 1;
}

sub get_database {
    my ($self) = @_;
    return $self->_param('database');
}

###########################################################################
#                       Data Specific Subroutines                         #
###########################################################################

sub get_all_transcripts {
    my ($self) = @_;
    
    my $query = 
        "SELECT t.feature_id, t.uniquename ".
        "FROM feature t, cvterm c ".
        "WHERE t.type_id = c.cvterm_id ".
        "AND c.name = 'transcript' ".
        "AND t.is_obsolete = 0";
    my $results = $self->do_select_query( $query );

    return $results;
}

sub get_all_rnas {
    my ($self) = @_;

    my $rna_type_string;
    map { $rna_type_string .= "'".$_."'," } @rna_types;
    $rna_type_string =~ s/,$//;

    my $query = 
        "SELECT f.feature_id, f.uniquename ".
        "FROM feature f, cvterm c ".
        "WHERE c.name IN ( $rna_type_string ) ".
        "AND c.cvterm_id = f.type_id ".
        "AND f.is_obsolete = 0";

    my $results = $self->do_select_query( $query );
    return $results;
}

sub get_cvterm_id {
    my ($self, $name) = @_;
    my $query = 
        "SELECT cvterm_id FROM cvterm where name = '$name'";
    my $results = $self->do_select_query( $query );
    return $results->[0]->[0] || undef;
}

## will find link between features in feature_relationship based on
## this data model:
# subject     rel_type     object
# transcript  derives_from gene
# polypeptide part_of      transcript
# CDS         derives_from transcript
# exon        part_of      transcript
# 
sub find_related {
    my ($self, $orig_feature_id, $class) = @_;
    my $retval;

    my $relationships = {
        'gene'        => [ 'subject', 'derives_from' ],
        'polypeptide' => [ 'object', 'part_of' ],
        'CDS'         => [ 'object', 'derives_from' ],
        'exon'        => [ 'object', 'part_of' ]
    };

    my $orig_feature_type = $self->get_feature_type( $orig_feature_id );
    
    unless( $orig_feature_type eq $class ) {
        
        my ($transrel, $reltype);
        if( $orig_feature_type eq 'transcript' ) {
            return undef unless( exists($relationships->{$class}) );
            ($transrel, $reltype) = @{$relationships->{$class}};
        } elsif( $class eq 'transcript' ) {
            return undef unless( exists($relationships->{$orig_feature_type}) );
            ($transrel, $reltype) = @{$relationships->{$orig_feature_type}};
        }
        
        if( defined( $transrel ) && defined( $reltype ) ) {

            my $query = 
                "SELECT f.feature_id ".
                "FROM feature f, cvterm c, cvterm r, feature_relationship fr ".
                "WHERE f.type_id = c.cvterm_id ".
                "AND c.name = '$class' ".
                "AND fr.type_id = r.cvterm_id ".
                "AND r.name = '$reltype' ";
            
            if( ( $transrel eq 'subject' && $orig_feature_type eq 'transcript' ) ||
                ( $transrel eq 'object'  && $class eq 'transcript' ) ) {
                $query .= 
                    "AND fr.subject_id = $orig_feature_id ".
                    "AND fr.object_id = f.feature_id";
            } elsif( ( $transrel eq 'object'  && $orig_feature_type eq 'transcript' ) ||
                     ( $transrel eq 'subject' && $class eq 'transcript' ) ) {
                $query .=
                    "AND fr.object_id = $orig_feature_id ".
                    "AND fr.subject_id = f.feature_id";
            }

            my $results = $self->do_select_query( $query );
            $retval = $results->[0]->[0] || undef;
        } else {

            #get the transcript id
            my $transcript_id = $self->find_related( $orig_feature_id, 'transcript' );

            if( !defined( $transcript_id ) ) {
                #maybe this is an rna (if it is, there's only one relationship
                my $query = 
                    "SELECT f.feature_id ".
                    "FROM feature f, cvterm c, cvterm r, feature_relationship fr ".
                    "WHERE c.cvterm_id = f.type_id ".
                    "AND c.name = '$class' ".
                    "AND r.name = 'derives_from' ".
                    "AND r.cvterm_id = fr.type_id ";

                if( grep ($orig_feature_type eq $_, @rna_types ) ) {
                    return undef unless( $class eq 'gene' );
                    $query .= 
                        "AND fr.subject_id = $orig_feature_id ".
                        "AND fr.object_id = f.feature_id";
                } elsif( $orig_feature_type eq 'gene' ) {
                    return undef unless( grep ( $class eq $_, @rna_types ) );
                    $query .=
                        "AND fr.subject_id = f.feature_id ".
                        "AND fr.object_id = $orig_feature_id";
                }
                my $results = $self->do_select_query( $query );
                $retval = $results->[0]->[0] || undef;
            } else {
            
                #get the related id 
                $retval = $self->find_related( $transcript_id, $class );
            }
        }

    } else {
        $retval = $orig_feature_id;
    }
    
    return $retval;
}

sub get_feature_type {
    my ($self, $feature_id) = @_;

    my $query = 
        "SELECT c.name FROM feature f, cvterm c ".
        "WHERE c.cvterm_id = f.type_id AND f.feature_id = $feature_id";
    my $results = $self->do_select_query( $query );

    if( @{$results} == 0 || @{$results->[0]} == 0 ) {
        die("Could not get feature type for feature $feature_id");
    }

    return $results->[0]->[0];
}
sub get_featureprop {
    my ($self, $feature_id, $prop) = @_;
    my $query = 
        "SELECT fp.value FROM featureprop fp, cvterm c".
        " WHERE fp.feature_id = $feature_id AND c.cvterm_id = fp.type_id ".
        "AND c.name = '$prop'";
    my $results = $self->do_select_query( $query );
    return $results->[0]->[0] || "";
}
sub get_all_featureprops {
    my ($self, $feature_id) = @_;
    my $query = 
        "SELECT c.name, fp.value FROM featureprop fp, cvterm c ".
        "WHERE fp.feature_id = $feature_id AND c.cvterm_id = fp.type_id";
    my $results = $self->do_select_query( $query );
    return $results;
}
sub get_all_organismprops {
    my ($self, $organism_id) = @_;
    my $query = 
        "SELECT cfp.name, fp.value ".
        "FROM organismprop fp, cvterm cfp ".
        "WHERE fp.organism_id = $organism_id ".
        "AND fp.type_id = cfp.cvterm_id";
    my $results = $self->do_select_query($query);
    return $results;
}
sub get_feature_attribute {
    my ($self, $feature_id, $att) = @_;
    my $query = 
        "SELECT $att FROM feature WHERE feature_id = $feature_id";
    my $results = $self->do_select_query( $query );
    return $results->[0]->[0];
}
sub get_residues {
    my ($self, $feature_id) = @_;
    return $self->get_feature_attribute( $feature_id, 'residues');
}
sub get_seqlen {
    my ($self, $feature_id) = @_;
    return $self->get_feature_attribute( $feature_id, 'seqlen' );
}
sub get_locus_id {
    my ($self, $feature_id, $db_name) = @_;

    if( !defined( $db_name ) ) {
        die("db_name is a required parameter to get_locus_id");
    }

    my $gene_id = $self->find_related( $feature_id, 'gene' );
    return undef unless( $gene_id );
    
    my $query = 
        "SELECT d.accession ".
        "FROM dbxref d, feature_dbxref fd, db ".
        "WHERE fd.dbxref_id = d.dbxref_id ".
        "AND db.db_id = d.db_id ".
        "AND db.name = '$db_name' ".
        "AND fd.feature_id = $gene_id";
    my $results = $self->do_select_query( $query );
    
    my $retval;
    unless( @{$results} == 0 || @{$results->[0]} == 0 ) {
        $retval = $results->[0]->[0];
    }

    return $retval;
}
sub get_external_locus_id {
    my ($self, $feature_id, $db_name) = @_;

    if( !defined( $db_name ) ) {
        die("db_name is a required parameter to get_external_locus_id");
    }

    my $gene_id = $self->find_related( $feature_id, 'gene' );
    return undef unless( $gene_id );
    
    my $query = 
        "SELECT d.accession ".
        "FROM dbxref d, feature_dbxref fd, db ".
        "WHERE fd.dbxref_id = d.dbxref_id ".
        "AND db.db_id = d.db_id ".
        "AND db.name = '$db_name' ".
        "AND fd.feature_id = $gene_id ".
        "AND d.version = 'current'";
    my $results = $self->do_select_query( $query );
    
    my $retval;
    unless( @{$results} == 0 || @{$results->[0]} == 0 ) {
        $retval = $results->[0]->[0];
    }

    return $retval;
}

sub get_feature_id_by_uniquename {
    my ($self, $uniquename) = @_;
    my $results = $self->do_select_query("SELECT feature_id FROM feature WHERE uniquename = '$uniquename'");
    
    if( @{$results} == 0 || @{$results->[0]} == 0 ) {
        die("Could not get feature_id for feature $uniquename");
    }
    return $results->[0]->[0];
}

sub get_feature_id_by_locus {
    my ($self, $locus, $db, $class) = @_;

    die("please pass in a locus to &get_feature_id_by_locus")if( !defined( $locus ) );
    die("please pass in a db name to &get_feature_id_by_locus") if( !defined( $db ) );

    my $query = 
        "SELECT f.feature_id ".
        "FROM feature f, dbxref d, feature_dbxref fd, db ".
        "WHERE f.feature_id = fd.feature_id ".
        "AND d.dbxref_id = fd.dbxref_id ".
        "AND db.name = '$db' ".
        "AND db.db_id = d.db_id ".
        "AND d.accession = '$locus'";

    my $results = $self->do_select_query( $query );
    if( @{$results} == 0 || @{$results->[0]} == 0 ) {
        die("Could not get feature_id for locus: $locus with given db name: $db");
    }
    
    my $retval;
    if( $class ) {
        $retval = $self->find_related( $results->[0]->[0], $class );
    } else {
        $retval = $results->[0]->[0];
    }
    return $retval;
}
sub get_organism_info {
    my ($self) = @_;
    my $query = 
        "SELECT o.abbreviation, o.genus, o.species, o.common_name, o.organism_id ".
        "FROM organism o ".
        "WHERE o.common_name != 'not known'";
    my $results = $self->do_select_query( $query );
    return $results;
}
sub get_molecule_info {
    my ($self) = @_;
    my $mol_string = join(", ", map { "'".$_."'" } @molecule_types );
    my $query = 
        "SELECT m.feature_id, m.uniquename, m.seqlen ".
        "FROM feature m, cvterm cm ".
        "WHERE m.type_id = cm.cvterm_id ".
        "AND cm.name IN ( $mol_string )";
    my $results = $self->do_select_query( $query );
    return $results;
}

sub get_molecule_and_coords {
    my ($self, $feature_id) = @_;
    my $query = 
        "SELECT m.uniquename, fl.fmin, fl.fmax, fl.strand ".
        "FROM feature m, featureloc fl ".
        "WHERE m.feature_id = fl.srcfeature_id ".
        "AND fl.feature_id = $feature_id";
    my $results = $self->do_select_query( $query );

    if( @{$results} == 0 || @{$results->[0]} == 0 ) {
        die("Could not get coordinates for feature $feature_id");
    }
    my $strand = $results->[0]->[3];
    if( $results->[0]->[3] == -1 ) {
        ($results->[0]->[1], $results->[0]->[2]) = ( $results->[0]->[2], $results->[0]->[1] );
    }
    return @{$results->[0]};
}
sub get_annotation {
    my ($self, $feature_id) = @_;
    my $transcript = $self->find_related( $feature_id, 'transcript' );

	my $transcript_uniquename = $self->get_uniquename_by_feature_id( $feature_id );
    
    #get the common_name
    my $common_name = $self->get_featureprop( $transcript, 'gene_product_name' );
    
    #get the gene symbol
    my $gene_sym = $self->get_featureprop( $transcript, 'gene_symbol' );

    #get the ec number
    my @ec_numbers = $self->get_ec_numbers( $transcript );

    #get the go ids
    my @go_ids  = $self->get_go_ids( $transcript );
    
    #get the tigr roles
    my @tigr_roles  = $self->get_tigr_roles( $transcript );
	
	# get product_name_source
    my $pn_source = $self->get_featureprop( $transcript, 'gene_product_name_source' );

	# get coord info
	my @coords = $self->get_coords_by_feature_id( $transcript );

    return ($common_name, $gene_sym, join(",",@ec_numbers), join(",", @go_ids), join(",", @tigr_roles), $transcript_uniquename, $pn_source, $coords[0], $coords[1], $coords[2] );
}

sub get_feature_cvterm_accessions {
    my ($self, $feature_id, @cvs) = @_;
    my $cv_string = join(", ", map { "'$_'" } @cvs );
    my $query = 
        "SELECT d.accession ".
        "FROM dbxref d, cv, feature_cvterm fc, cvterm c ".
        "WHERE cv.name IN ($cv_string) ".
        "AND cv.cv_id = c.cv_id ".
        "AND d.dbxref_id = c.dbxref_id ".
        "AND c.cvterm_id = fc.cvterm_id ".
        "AND fc.feature_id = $feature_id";
    my $results = $self->do_select_query( $query );

    my @retval = ();
    if( @{$results} != 0 ) {
        map { push(@retval, $_->[0]) } @{$results};
    }
    return @retval;
}

sub get_ec_numbers {
    my ($self, $feature_id) = @_;
    my @retval;
    my @results = $self->get_feature_cvterm_accessions( $feature_id, 'EC' );

    return () unless( @results > 0 );
    
    my $accession_string = join( ", ", map { "'$_'" } @results );
    my $query = 
        "SELECT d2.accession ".
        "FROM dbxref d1, dbxref d2, cvterm_dbxref cd, cvterm c, cv, db ".
        "WHERE cv.name = 'EC' ".
        "AND c.cv_id = cv.cv_id ".
        "AND c.cvterm_id = cd.cvterm_id ".
        "AND cd.dbxref_id = d2.dbxref_id ".
        "AND d1.accession IN ($accession_string) ".
        "AND d1.dbxref_id = c.dbxref_id ".
        "AND d2.db_id = db.db_id ".
        "AND db.name = 'EC'";
    my $results = $self->do_select_query( $query );

    @retval = map { $_->[0] } @{$results};
    return @retval;
}
sub get_go_ids {
    my ($self, $feature_id) = @_;
    return $self->get_feature_cvterm_accessions( $feature_id, ('GO','function','process','component','biological_process','molecular_function','cellular_component'));
}
sub get_tigr_roles {
    my ($self, $feature_id) = @_;
    my @retval;
    my @results = $self->get_feature_cvterm_accessions( $feature_id, 'TIGR_role' );

    return () unless( @results > 0 );
    
    my $accession_string = join( ", ", map { "'$_'" } @results );
    my $query = 
        "SELECT d2.accession ".
        "FROM dbxref d1, dbxref d2, cvterm_dbxref cd, cvterm c, cv, db ".
        "WHERE cv.name = 'TIGR_role' ".
        "AND c.cv_id = cv.cv_id ".
        "AND c.cvterm_id = cd.cvterm_id ".
        "AND cd.dbxref_id = d2.dbxref_id ".
        "AND d1.accession IN ($accession_string) ".
        "AND d1.dbxref_id = c.dbxref_id ".
        "AND d2.db_id = db.db_id ".
        "AND db.name = 'TIGR_role'";
    my $results = $self->do_select_query( $query );

    @retval = map { $_->[0] } @{$results};
    return @retval;
}

sub get_uniquename_by_feature_id {
    my ($self, $feature_id ) = @_;
    my $query = 
        "SELECT uniquename FROM feature WHERE feature_id = $feature_id";

    my $results = $self->do_select_query( $query );
    if( @{$results} == 0 || @{$results->[0]} == 0 ) {
        die("Could not get uniquename for feature with feature_id: $feature_id");
    }
    return $results->[0]->[0];
}

sub get_coords_by_feature_id {
    my ($self, $feature_id ) = @_;
	my @retval;
	
    my $query = 
        "SELECT fl.fmin, fl.fmax, fl.strand FROM featureloc fl WHERE fl.feature_id = $feature_id";

    my $results = $self->do_select_query( $query );
    if( @{$results} == 0 || @{$results->[0]} == 0 ) {
        die("Could not get coords for feature with feature_id: $feature_id");
    }
	
	@retval = map { $_->[0], $_->[1], $_->[2] } @{$results};
	return @retval;
}

sub feature_is_complete {
    my ($self, $feature_id) = @_;
    
    my $transcript_id = $self->find_related( $feature_id, 'transcript' );
    return unless( $transcript_id );
    my $query = 
        "SELECT fp.value FROM featureprop fp, cvterm c ".
        "WHERE c.name = 'completed_by' ".
        "AND c.cvterm_id = fp.type_id ".
        "AND fp.feature_id= ?";
    
    my $retval = 0;
    my $results = $self->do_select_query( $query, $feature_id );
    if( @{$results} > 0 ) {
        $retval = $results->[0]->[0];
    }
    $retval;
}

###########################################################################
#                       Private Subroutines                               #
###########################################################################

sub _connect {
    my ($self) = @_;
    my ($db, $host,$user, $pass) = ($self->_param('database'), $self->_param('host'), 
                                    $self->_param('username'), $self->_param('password') );
    my $retval = DBI->connect("dbi:mysql:host=$host;packetSize=8092", $user, $pass,
                              { RaiseError => 1 }
                              ) or
        die("Could not connect to server:".DBI->errstr);

    $self->_param('dbh', $retval);
    if( defined( $db ) ) {
        $self->change_db( $db );
    }
}

sub _init {
    my ($self, %args) = @_;

    my @reqs = qw(username password host);
    foreach my $req ( @reqs ) {
        if( exists( $args{$req} ) ) {
            $self->_param( $req, $args{$req} );
        } else {
            die("$req is a required parameter to ".__PACKAGE__."::new");
        }
    }

    my @opts = qw(database);
    foreach my $opt ( @opts ) {
        $self->_param( $opt, $args{$opt} ) if( exists( $args{$opt} ) );
    }

    $self->_connect();
}

sub _param {
    my ($self, $param, $value) = @_;
    if( defined( $value ) ) {
        $self->{$param} = $value;
    }
    return $self->{$param};
}

1==1;
