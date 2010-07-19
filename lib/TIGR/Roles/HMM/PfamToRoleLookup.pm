package TIGR::Roles::HMM::PfamToRoleLookup;

use strict;
use warnings;
use MLDBM 'DB_File';
use Fcntl qw( O_RDONLY );
use Data::Dumper;

my $default_pfam2role_dbfile = "/pfam2role.db";
my $default_pfam2role_tab = "/hmm/PFAM_role_id.txt";

sub new {
    my ($class, %args) = @_;
    my $self = {};
    bless($self,$class);
    $self->_init(%args);
    return $self;
}

sub DESTROY {
    my ($self) = @_;
    my $pfam_lookup = $self->tied_lookup_hash();
    untie %{$pfam_lookup};
}

sub pfam2tigr_role {
    my ($self, $pfam_acc) = @_;
    my $lookup = $self->tied_lookup_hash();
    return $lookup->{$pfam_acc} || [];
}

sub _init {
    my ($self, %args) = @_;

    if( $args{'roles_db_dir'} ) {
        $self->roles_db_dir($args{'roles_db_dir'});
    }

    if( $args{'pfam2role_dbfile'} ) {
        $self->pfam2role_dbfile( $args{'pfam2role_dbfile'} );
    } else {
        $self->pfam2role_dbfile( $self->roles_db_dir().$default_pfam2role_dbfile )
            if( $self->roles_db_dir() );
    }

    if( $args{'pfam2role_tab'} ) {
        $self->pfam2role_tab( $args{'pfam2role_tab'} );
    } else {
        $self->pfam2role_tab( $self->roles_db_dir().$default_pfam2role_tab )
            if( $self->roles_db_dir() );
    }

    #die if we don't have the files we need
    die("Please provide paths to required files by using roles_db_dir, pfam2role_dbfile and/or pfam2role_tab options")
        unless( $self->pfam2role_tab && $self->pfam2role_dbfile );

    #Tie a hash
    my %pfam2role_lookup;
    tie( %pfam2role_lookup, 'MLDBM', $self->pfam2role_dbfile, O_RDONLY );
    $self->tied_lookup_hash( \%pfam2role_lookup );

    if( ! -e $self->pfam2role_dbfile || $args{'force'} ) {
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

    my $pfam_txt = $self->pfam2role_tab();
    open( IN, "< $pfam_txt") or die("Could not open $pfam_txt ($!)");

    my $lookup = $self->tied_lookup_hash();
    
    #Reset the hash
    %{$lookup} = ();

    my $tmp;
    while(<IN>) {
        my ($pfam,$role_id) = split(/\s+/);
        push(@{$tmp->{$pfam}}, $role_id );
    }

    %{$lookup} = %{$tmp};
    
    close(IN);
}

sub roles_db_dir {
    my ($self, $value) = @_;
    return $self->_param('roles_db_dir', $value);
}

sub pfam2role_dbfile {
    my ($self, $value) = @_;
    return $self->_param('pfam2role_dbfile', $value);
}

sub pfam2role_tab {
    my ($self, $value) = @_;
    return $self->_param('pfam2role_tab', $value);
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
        return $self->{"_$name"} || undef;
    }
}
1;
