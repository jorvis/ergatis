package TIGR::Roles::Omnium::OmniumToRoleLookup;

use strict;
use warnings;
use MLDBM 'DB_File';
use Carp;
use Data::Dumper;

my($LOCUS,$ROLE_ID) = (1..2);

sub new {
    my ($class, %args) = @_;
    my $self = {};
    bless($self, $class);
    $self->_init(%args);
    return $self;
}

sub DESTROY {
    my ($self) = @_;
    my $lookup = $self->tied_lookup_hash();
    untie %{$lookup};
}

sub _init {
    my ($self, %args) = @_;

    if( $args{'roles_db_dir'} ) {
        $self->roles_db_dir($args{'roles_db_dir'});
    } else {
        $self->roles_db_dir("/usr/local/projects/db/tigr_roles");
    }

    if( $args{'omnium2role_dbfile'} ) {
        $self->omnium2role_dbfile( $args{'omnium2role_dbfile'} );
    } else {
        $self->omnium2role_dbfile( $self->roles_db_dir()."/omnium2role.db" );
    }

    if( $args{'bcp_role_link'} ) {
        $self->bcp_role_link( $args{'bcp_role_link'} );
    } else {
        $self->bcp_role_link( $self->roles_db_dir()."/omnium/bcp/bcp_role_link" );
    }

    #Tie a hash
    my %omnium2role_lookup;
    tie( %omnium2role_lookup, 'MLDBM', $self->omnium2role_dbfile );
    $self->tied_lookup_hash( \%omnium2role_lookup );

    if( ! -e $self->omnium2role_dbfile || $args{'force'} ) {
        $self->_create_lookup();
    }

}

sub get_count {
    my ($self) = @_;
    my $lookup = $self->tied_lookup_hash();
    return scalar(keys(%{$lookup}));
}

sub tied_lookup_hash {
    my ($self, $value) = @_;
    return $self->_param('tied_lookup_hash', $value);
}

sub roles_db_dir {
    my ($self, $value) = @_;
    return $self->_param('roles_db_dir', $value );
}
sub omnium2role_dbfile {
    my ($self, $value) = @_;
    return $self->_param('omnium2role_dbfile', $value);
}
sub bcp_role_link {
    my ($self, $value) = @_;
    return $self->_param('bcp_role_link', $value);
}

sub _create_lookup {
    my ($self) = @_;

    my $bcp_role_link = $self->bcp_role_link();

    open( IN, "< $bcp_role_link") or croak("Can't open $bcp_role_link $!");

    my $lookup = $self->tied_lookup_hash;

    #Reset the tied hash
    %{$lookup} = ();

    my $tmp;
    while(<IN>) {
        my @cols = split(/\s+/);
        push( @{ $tmp->{ $cols[$LOCUS] } }, $cols[$ROLE_ID] );
    }

    %{$lookup} = %{$tmp};

    close(IN);
}

sub omnium2tigr_role {
    my ($self, $omnium_locus) = @_;
    my $lookup = $self->tied_lookup_hash();
    return $lookup->{$omnium_locus} || [];
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
