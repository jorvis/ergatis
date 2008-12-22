package TIGR::Roles::NameLookup;

use strict;
use warnings;
use MLDBM "DB_File";
use Data::Dumper;

my $default_roles_db_dir = "/usr/local/projects/db/tigr_roles";
my $default_role2name_dbfile = "/role_data.db";
my $default_role2name_bcp = "/omnium/bcp/bcp_egad_roles";

my ($ROLE_ID, $COMPARTMENT, $MAINROLE, $SUPERROLE, $SUB1ROLE, $SUB2ROLE, $SUB3ROLE, $SUB4ROLE) = (0..7);

sub new {
    my ($class, %args) = @_;
    my $self = {};
    bless($self,$class);
    $self->_init(%args);
    return $self;
}

sub DESTROY {
    my ($self) = @_;
    my $roles_data_lookup = $self->tied_lookup_hash();
    untie %{$roles_data_lookup};
}

sub _init {
    my ($self, %args) = @_;

    if( $args{'roles_db_dir'} ) {
        $self->roles_db_dir($args{'roles_db_dir'});
    } else {
        $self->roles_db_dir($default_roles_db_dir);
    }

    if( $args{'role2name_dbfile'} ) {
        $self->role2name_dbfile( $args{'role2name_dbfile'} );
    } else {
        $self->role2name_dbfile( $self->roles_db_dir().$default_role2name_dbfile );
    }

    if( $args{'role2name_bcp'} ) {
        $self->role2name_bcp( $args{'role2name_bcp'} );
    } else {
        $self->role2name_bcp( $self->roles_db_dir().$default_role2name_bcp );
    }

    #Tie a hash
    my %role_data_lookup;
    tie( %role_data_lookup, 'MLDBM', $self->role2name_dbfile );
    $self->tied_lookup_hash( \%role_data_lookup );

    if( ! -e $self->role2name_dbfile || $args{'force'} ) {
        $self->_create_lookup();
    }

}

sub get_count {
    my ($self) = @_;
    my $lookup = $self->tied_lookup_hash;
    return scalar( keys( %{$lookup} ) );
}

sub _create_lookup {
    my ($self) = @_;

    my $egad_roles_bcp = $self->role2name_bcp();
    open( IN, "< $egad_roles_bcp") or die("Could not open $egad_roles_bcp ($!)");

    my $lookup = $self->tied_lookup_hash();
    
    #Reset the hash
    %{$lookup} = ();

    my $tmp;
    while(<IN>) {
        chomp;
        my @cols = split(/\t/);

        my $role_id = $cols[$ROLE_ID];
        my %this_role = ( 'compartment' => $cols[$COMPARTMENT],
                          'main_role'   => $cols[$MAINROLE],
                          'super_role'  => $cols[$SUPERROLE],
                          'sub1role'    => $cols[$SUB1ROLE],
                          'sub2role'    => $cols[$SUB2ROLE],
                          'sub3role'    => $cols[$SUB3ROLE],
                          'sub4role'    => $cols[$SUB4ROLE], );
        
        $tmp->{$role_id} = \%this_role;   
    }
    
    %{$lookup} = %{$tmp};
    
    close(IN);
}

## Lookup subroutines
sub compartment {
    my ($self, $role_id) = @_;
    return $self->_lookup( $role_id, 'compartment' );
}
sub main_role {
    my ($self, $role_id) = @_;
    return $self->_lookup( $role_id, 'main_role' );
}
sub super_role {
    my ($self, $role_id) = @_;
    return $self->_lookup( $role_id, 'super_role' );
}
sub sub1role {
    my ($self, $role_id) = @_;
    return $self->_lookup( $role_id, 'sub1role' );
}
sub sub2role {
    my ($self, $role_id) = @_;
    return $self->_lookup( $role_id, 'sub2role' );
}
sub sub3role {
    my ($self, $role_id) = @_;
    return $self->_lookup( $role_id, 'sub3role' );
}
sub sub4role {
    my ($self, $role_id) = @_;
    return $self->_lookup( $role_id, 'sub4role' );
}
sub _lookup {
    my ($self, $role_id, $data_class) = @_;
    my $lookup = $self->tied_lookup_hash();
    my $retval = undef;

    if( defined($lookup->{$role_id}) ) {
        $retval = $lookup->{$role_id}->{$data_class};
    }
    
    return $retval;
}

## Object accessor
sub roles_db_dir {
    my ($self, $value) = @_;
    return $self->_param('roles_db_dir', $value);
}

sub role2name_dbfile {
    my ($self, $value) = @_;
    return $self->_param('role2name_dbfile', $value);
}

sub role2name_bcp {
    my ($self, $value) = @_;
    return $self->_param('role2name_bcp', $value);
}

sub tied_lookup_hash {
    my ($self, $value) = @_;
    return $self->_param('tied_lookup_hash', $value);
}

sub _param {
    my ($self,$name,$value) = @_;
    if( $value ) {
        $self->{"_$name"} = $value;
    } else {
        return $self->{"_$name"};
    }
}
1;
