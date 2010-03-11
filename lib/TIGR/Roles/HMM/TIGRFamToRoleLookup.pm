package TIGR::Roles::HMM::TIGRFamToRoleLookup;

use strict;
use warnings;
use MLDBM 'DB_File';
use Fcntl qw( O_RDONLY );
use Data::Dumper;

my $default_roles_db_dir = "/usr/local/projects/db";
my $default_tigrfam2role_dbfile = "/tigr_roles/tigrfam2role.db";
my $default_tigrfam2role_tab = "/TIGRFAMs/TIGRFAMS_ROLE_LINK";

sub new {
    my ($class, %args) = @_;
    my $self = {};
    bless($self,$class);
    $self->_init(%args);
    return $self;
}

sub tigrfam2tigr_role {
    my ($self, $tigrfam_acc ) = @_;
    my $lookup = $self->tied_lookup_hash();
    return $lookup->{$tigrfam_acc} || [];
}

sub DESTROY {
    my ($self) = @_;
    my $tigrfam_lookup = $self->tied_lookup_hash();
    untie %{$tigrfam_lookup};
}

sub _init {
    my ($self, %args) = @_;

    if( $args{'roles_db_dir'} ) {
        $self->roles_db_dir($args{'roles_db_dir'});
    } else {
        $self->roles_db_dir($default_roles_db_dir);
    }

    if( $args{'pfam2role_dbfile'} ) {
        $self->tigrfam2role_dbfile( $args{'tigrfam2role_dbfile'} );
    } else {
        $self->tigrfam2role_dbfile( $self->roles_db_dir().$default_tigrfam2role_dbfile );
    }

    if( $args{'tigrfam2role_tab'} ) {
        $self->tigrfam2role_tab( $args{'tigrfam2role_tab'} );
    } else {
        $self->tigrfam2role_tab( $self->roles_db_dir().$default_tigrfam2role_tab );
    }

    #Tie a hash
    my %tigrfam2role_lookup;
    tie( %tigrfam2role_lookup, 'MLDBM', $self->tigrfam2role_dbfile, O_RDONLY );
    $self->tied_lookup_hash( \%tigrfam2role_lookup );

    if( ! -e $self->tigrfam2role_dbfile || $args{'force'} ) {
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

    my $tigrfam_txt = $self->tigrfam2role_tab();
    open( IN, "< $tigrfam_txt") or die("Could not open $tigrfam_txt ($!)");

    my $lookup = $self->tied_lookup_hash();
    
    #Reset the hash
    %{$lookup} = ();

    my $tmp;
    while(<IN>) {
        my ($tigrfam,$role_id) = split(/\s+/);
        push( @{$tmp->{$tigrfam}}, $role_id );
    }

    %{$lookup} = %{$tmp};

    close(IN);
}

sub roles_db_dir {
    my ($self, $value) = @_;
    return $self->_param('roles_db_dir', $value);
}

sub tigrfam2role_dbfile {
    my ($self, $value) = @_;
    return $self->_param('tigrfam2role_dbfile', $value);
}

sub tigrfam2role_tab {
    my ($self, $value) = @_;
    return $self->_param('tigrfam2role_tab', $value);
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


