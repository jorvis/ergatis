package TIGR::Roles::HMM::TIGRFamToRoleLookup;

use strict;
use warnings;
use MLDBM 'DB_File';
use Fcntl qw( O_RDONLY );
use Data::Dumper;

my $default_tigrfam2role_tab = "/TIGRFAMS_ROLE_LINK";

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

    if( $args{'tigrfams_dir'} ) {
        $self->tigrfams_dir( $args{'tigrfams_dir'} );
    }

    if( $args{'tigrfam2role_tab'} ) {
        $self->tigrfam2role_tab( $args{'tigrfam2role_tab'} );
    } else {
        $self->tigrfam2role_tab( $self->tigrfams_dir().$default_tigrfam2role_tab );
    }

	# Create lookup hash
    $self->_create_lookup();


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
    
    #Reset the hash
    my $lookup;    
    %{$lookup} = ();

    my $tmp;
    while(<IN>) {
        my ($tigrfam,$role_id) = split(/\s+/);
        push( @{$tmp->{$tigrfam}}, $role_id );
    }

    %{$lookup} = %{$tmp};

	#set the lookup hash 
	$self->tied_lookup_hash($lookup);
	
    close(IN);
}

sub tigrfams_dir {
    my ($self, $value) = @_;
    return $self->_param('tigrfams_dir', $value);
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
        return $self->{"_$name"} || undef;
    }
}
1;


