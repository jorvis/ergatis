package IPD::IPDObject::Project;


=head1 NAME 

IPD::IPDObject::Project - Module for managing attributes of Project-type IPD objects

=cut

use IPD::Client;
use strict;
use warnings;
use base qw(IPD::IPDObject);

our $AUTOLOAD;

my $extension = "";

my $attributes = {
	'name' => -1, 
	'contact' => -1,
	'description' => -1, 
	'id' => -1,
	'lims_project_id' => -1,
	'status' => -1,
	'is_active' => -1,
	'update_at' => -1,
	'created_at' => -1, 
	'project_properties' => -1,
	'fundings' => -1
};


sub new {
    my ($class, %args) = @_;
    my $self = $class->SUPER::new( %args );
    bless($self, $class);
    return $self;
}

sub AUTOLOAD {
    my ($self, $value) = @_;
    my ($oper, $attr) = ($AUTOLOAD =~ /(get_|set_)(\w+)$/);
    return if( $AUTOLOAD =~ /DESTROY/ );

    unless ($oper && $attr) {
	die "Method name $AUTOLOAD is not in the recognized form (get|set)_attribute\n";
    }
    unless (exists ${$attributes}{$attr}) {
	die "No such attribute '$attr' exists for the class " . ref($self) . "\n";
    }

    no strict 'refs';
    # adds $AUTOLOAD name to symbol table for quicker repeated calls to same function
    if ($oper eq 'get_') {	
#	*{$AUTOLOAD} = sub { shift->{$attr} };
    } elsif ($oper eq 'set_') {
#	*{$AUTOLOAD} = sub { shift->{$attr} = shift };
        $self->{$attr} = $value;	# if we want to 'set' the new value
    }
    use strict 'refs';

    return $self->{$attr};
}

### Attribute Get Methods that cannot be easily accessed with AUTOLOAD ###
sub get_funding_id {
    my $self = shift;
    return $self->{'fundings'}->{'funding'}->{'id'};
}

sub get_funding_name {
    my $self = shift;
    return $self->{'fundings'}->{'funding'}->{'name'};
}

sub get_funding_start_date {
    my $self = shift;
    return $self->{'fundings'}->{'funding'}->{'start_date'};
}

sub get_funding_end_date {
    my $self = shift;
    return $self->{'fundings'}->{'funding'}->{'end_date'};
}

sub get_funding_agent {
    my $self = shift;
    return $self->{'fundings'}->{'funding'}->{'funding_agent'};
}

sub get_funding_status {
    my $self = shift;
    return $self->{'fundings'}->{'funding'}->{'status'};
}

sub get_funding_description {
    my $self = shift;
    return $self->{'fundings'}->{'funding'}->{'description'};
}

sub get_funding_created_at {
    my $self = shift;
    return $self->{'fundings'}->{'funding'}->{'created_at'};
}

sub get_funding_updated_at {
    my $self = shift;
    return $self->{'fundings'}->{'funding'}->{'updated_at'};
}

sub get_funding_first_name {
    my $self = shift;
    return $self->{'fundings'}->{'funding'}->{'contact'}->{'first_name'};
}

sub get_funding_last_name {
    my $self = shift;
    return $self->{'fundings'}->{'funding'}->{'contact'}->{'last_name'};
}

sub get_funding_email {
    my $self = shift;
    return $self->{'fundings'}->{'funding'}->{'contact'}->{'email'};
}

sub get_funding_username {
    my $self = shift;
    return $self->{'fundings'}->{'funding'}->{'contact'}->{'username'};
}

sub get_property_id {
    my $self = shift;
    return $self->{'project_properties'}->{'project_property'}->{'id'};
}

sub get_property_name {
    my $self = shift;
    return $self->{'project_properties'}->{'project_property'}->{'name'};
}

sub get_property_value {
    my $self = shift;
    return $self->{'project_properties'}->{'project_property'}->{'value'};
}

sub get_property_created_at {
    my $self = shift;
    return $self->{'project_properties'}->{'project_property'}->{'created_at'};
}

sub get_property_updated_at {
    my $self = shift;
    return $self->{'project_properties'}->{'project_property'}->{'updated_at'};
}

sub get_updated_at {	# Only projects XML has 'update_at' instead of 'updated_at' so this is to avoid confusion
    my $self = shift;
    return $self->{'update_at'};
}

#Methods we do not want to access with AUTOLOAD
sub set_id {
    my $self = shift;
    my $value = shift;
    die ("Not allowed to change the ID: $!\n") if (defined $self -> {'id'});
    $self -> {'id'} = $value;
}

#Class and instance methods using this class

=item IPD::IPDObject::Project->get_all_projects()

B<Description:> Retrieves all project information from IPD and stores into an array of objects
B<Parameters:> None
B<Returns:> An array reference of IPD::IPDObject::Project objects

=cut

sub get_all_projects {
    my $class = shift;
    die ("Must call 'get' methods in the form of ", ref($class), "->method() : $!\n") if ref $class;
    my @selves;
    $extension = "projects.xml";
    my $elem = $class->SUPER::XML2Object( IPD::Net->send_request( IPD::Client->get_client, $extension, "GET" ) );
    if (ref($elem->{'project'}) eq 'HASH') {
	push @selves, $class->new(%{$elem->{'project'}});
    } else {
    	foreach my $proj ( @{$elem->{'project'}} ) {
	    push @selves,$class->new(%{$proj}); 
        }
    }
    return \@selves;
}

=item IPD::IPDObject::Project->get_project($id)

B<Description:> Retrieves project information of a single project ID from IPD and stores into an object
B<Parameters:> The ID of the project to be queried.
B<Returns:> An IPD::IPDObject::Project object

=cut
sub get_project {
    my $class = shift;
    die ("Must call 'get' methods in the form of ", ref($class), "->method() : $!\n") if ref $class;
    my $id = shift;
    $extension ="projects/$id.xml";
    my $elem = $class->SUPER::XML2Object( IPD::Net->send_request( IPD::Client->get_client,, $extension, "GET" ) );
    my $proj = $class->new(%{$elem});
    return $proj;
}

1;
