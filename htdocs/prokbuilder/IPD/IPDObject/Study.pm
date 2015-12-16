package IPD::IPDObject::Study;

=head1 NAME 

IPD::IPDObject::Study - Module for managing attributes of Study-type IPD objects

=cut

use IPD::Client;
use IPD::IPDObject::Sample;
use IPD::IPDObject::StudyProperty;
use IPD::IPDObject::StudyStage;

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
	'project_id' => -1,
	'status' => -1,
	'type' => -1,
	'updated_at' => -1,
	'created_at' => -1,
	'genome_project_id' => -1,
	'targeted_finishing_goal' => -1, 
	'genome_project_id' => -1,

	'sample' => -1,
	'study_stage' => -1,
	'study_stages' => -1,
	'study_property' => -1,
	'study_properties' => -1
};

sub new {
    my ($class, %args) = @_;
    my $self = $class->SUPER::new(%args);
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

#Methods we do not want to access with AUTOLOAD
sub set_id {
    my $self = shift;
    my $value = shift;
    die ("Not allowed to change the ID: $!\n") if (defined $self -> {'id'});
    $self -> {'id'} = $value;
}

sub set_type {
    my $self = shift;
    my $value = shift;
    die ("Not allowed to change the type: $!\n") if (defined $self->{'type'});
    $self->{'type'} = $value;
}

#Manipulating other objects through this object
sub get_study_stages {
    my $self = shift;
    return IPD::IPDObject::StudyStage->get_study_stages_by_study( $self->{'id'} );
}

sub get_sample {
    my $self = shift;
    return IPD::IPDObject::Sample->get_sample_by_study( $self->{'id'} );
}

sub get_study_properties {
    my $self = shift;
    return IPD::IPDObject::StudyProperty->get_study_properties_by_study( $self->{'id'} );
}

#Class and instance methods using this class

=item IPD::IPDObject::Study->get_study($id)

B<Description:> Retrieves study information of a single study ID from IPD and stores into an object
B<Parameters:> The ID of the study to be queried.
B<Returns:> An IPD::IPDObject::Study object

=cut

sub get_study {
    my $class = shift;
    die ("Must call 'get' methods in the form of ", ref($class), "->method() when returning an object of this type: $!\n") if ref $class;
    my $id = shift;

    $extension ="studies/$id.xml";
    my $elem = $class->SUPER::XML2Object(IPD::Net->send_request(IPD::Client->get_client, $extension, "GET"));
    my $study = $class->new(%{$elem});
    return $study;
}


=item IPD::IPDObject::Study->get_all_studies

B<Description:> Retrieves all study information from IPD and stores into an array of objects
B<Parameters:> None.
B<Returns:> An array of IPD::IPDObject::Study objects
B<Note:>  Because of the sheer amount of studies currently present in IPD, this command may give a "time out" error due to taking too long (the time limit set is 999 seconds).  If this error occurs, contact me (sadkins@som.umaryland.edu) and I will make the limit longer

=cut

sub get_all_studies {
    my $class = shift;
    die ("Must call 'get' methods in the form of ", ref($class), "->method() when returning an object of this type: $!\n") if ref $class;
    my @selves;

    $extension = "studies.xml";
    my $elem = $class->SUPER::XML2Object(IPD::Net->send_request(IPD::Client->get_client, $extension, "GET"));
    if (ref($elem->{'study'}) eq 'HASH') {
	push @selves, $class->new(%{$elem->{'study'}});
    } else {
    	foreach my $study ( @{$elem->{'study'}} ) {
	    push @selves, $class->new(%{$study}); 
        }
    }
    return \@selves;
}

=item IPD::IPDObject::Study>get_all_studies_by_project($id)

B<Description:> Retrieves all study information associated with a single project ID from IPD and stores into an array of objects
B<Parameters:> The ID of the project to be queried.
B<Returns:> An array of IPD::IPDObject::Study objects

=cut

sub get_all_studies_by_project {
    my $class = shift;
    die ("Must call 'get' methods in the form of ", ref($class), "->method() when returning an object of this type: $!\n") if ref $class;
    my $id = shift;
    my @selves;

    $extension = "studies.xml?project_id=$id";
    my $elem = $class->SUPER::XML2Object(IPD::Net->send_request(IPD::Client->get_client, $extension, "GET"));
    if (ref($elem->{'study'}) eq 'HASH') {
	push @selves, $class->new(%{$elem->{'study'}});
    } else {
    	foreach my $study ( @{$elem->{'study'}} ) {
	    push @selves, $class->new(%{$study}); 
        }
    }
    return \@selves;
}

=item IPD::IPDObject::Study->get_all_studies_by_lims_id($id)

B<Description:> Retrieves all study information associated with a single LIMS ID from IPD and stores into an array of objects
B<Parameters:> The LIMS ID to be queried.
B<Returns:> An array of IPD::IPDObject::Study objects

=cut

sub get_all_studies_by_lims_id {
    my $class = shift;
    die ("Must call 'get' methods in the form of ", ref($class), "->method() when returning an object of this type: $!\n") if ref $class;
    my $id = shift;
    my @selves;

    $extension = "studies.xml?lims_project_id=$id";
    my $elem = $class->SUPER::XML2Object(IPD::Net->send_request(IPD::Client->get_client, $extension, "GET"));
    if (ref($elem->{'study'}) eq 'HASH') {
	push @selves, $class->new(%{$elem->{'study'}});
    } else {
    	foreach my $study ( @{$elem->{'study'}} ) {
	    push @selves, $class->new(%{$study}); 
        }
    }
    return \@selves;
}

=item IPD::IPDObject::Study->get_studies_by_lims_master($id)

B<Description:> Retrieves all studies associated with a single LIMS Master ID from IPD and stores into an object
B<Parameters:> The LIMS Master ID to be queried.
B<Returns:> An array of IPD::IPDObject::Study objects

=cut

sub get_studies_by_lims_master {
    my $class = shift;
    die ("Must call 'get' methods in the form of ", ref($class), "->method() when returning an object of this type: $!\n") if ref $class;
    my $id = shift;
    my @selves;

    $extension = "studies.xml?lims_sample_id=$id";
    my $elem = $class->SUPER::XML2Object(IPD::Net->send_request(IPD::Client->get_client, $extension, "GET"));
    if (ref($elem->{'study'}) eq 'HASH') {
	push @selves, $class->new(%{$elem->{'study'}});
    } else {
    	foreach my $study ( @{$elem->{'study'}} ) {
	    push @selves, $class->new(%{$study}); 
        }
    }
    return \@selves;
}

=item $obj->create_study()

B<Description:> Creates a new Study entry in IPD based on the provided Study object
B<Parameters:> An IPD::IPDObject::Study object
B<Returns:> The same IPD::IPDObject::Study object

=cut

sub create_study {
    my $self = shift;
    die ("Must call 'create|save|delete' methods in the form of ", '$object->method()', " : $!\n") unless ref $self;

    foreach my $name ('name', 'project_id', 'type', 'status') {
	if (!exists(${$self}{$name})) {
	    die "ERROR:  name, project_id, type, and status elements are required. $!\n";	
	}
    }
    $extension = "studies.xml";
    my $xml =  ref($self)->SUPER::construct_xml($self, "study");
    my $temp = ref($self)->SUPER::XML2Object(IPD::Net->send_request(IPD::Client->get_client, $extension, "POST", $xml));
    my $id = $temp->{'id'};
    $extension = "studies/$id.xml";
    my $elem = ref($self)->SUPER::XML2Object(IPD::Net->send_request(IPD::Client->get_client, $extension, "GET"));
    $self = ref($self)->new(%{$elem});
    return $self;
}

=item $obj->save_study()

B<Description:> Updates a Study entry in IPD based on the provided Study object.  The object must have an 'id' attribute associated with an existing study in IPD
B<Parameters:> An IPD::IPDObject::Study object
B<Returns:> The same IPD::IPDObject::Study object

=cut

sub save_study {
    my $self = shift;
    die ("Must call 'create|save|delete' methods in the form of ", '$object->method()', " : $!\n") unless ref $self;
    my $id;

    exists $self->{'id'} ? $id = $self->{'id'} : die "No Study ID present.  Perhaps you wanted to use 'create_study' method? $!\n";
    $extension = "studies/$id.xml";
    my $xml =  ref($self)->SUPER::construct_xml($self, "study");
    IPD::Net->send_request(IPD::Client->get_client, $extension, "PUT", $xml);
    my $elem = ref($self)->SUPER::XML2Object(IPD::Net->send_request(IPD::Client->get_client, $extension, "GET"));
    $self = ref($self)->new(%{$elem});
    return $self;
}


=item $obj->delete_study()

B<Description:> Deletes a Study entry from IPD, based on the ID attribute found in the Study object.  Also clears all attributes from the object.
B<Parameters:> An IPD::IPDObject::Study object
B<Returns:> Nothing

=cut

sub delete_study {
    my $self = shift;
    die ("Must call 'create|save|delete' methods in the form of ", '$object->method()', " : $!\n") unless ref $self;

    die "No Study ID attribute present for delete command: $!\n" if (! exists $self->{'id'});
    my $id = $self->{'id'};
    $extension = "studies/$id.xml";
    IPD::Net->send_request(IPD::Client->get_client, $extension, "DELETE");
    $self = ();
    return;
}

1;
