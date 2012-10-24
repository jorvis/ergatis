package IPD::IPDObject::StudyStage;

=head1 NAME 

IPD::IPDObject::StudyStage - Module for managing attributes of StudyStage-type IPD objects

=cut

use IPD::Client;
use IPD::IPDObject::StudyStageProperty;

use strict;
use warnings;
use base qw(IPD::IPDObject);

our $AUTOLOAD;

my $extension = "";

my $attributes = {
	'contact' => -1,
	'description' => -1, 
	'id' => -1,
	'study_id' => -1,
	'status' => -1,
	'type' => -1,
	'updated_at' => -1,
	'created_at' => -1,
	'start_date' => -1,
	'target_date' => -1,
	'end_date' => -1,
	'is_active' => -1,
	'lims_work_order_id' => -1,

	'study_stage_properties' => -1,
	'study_stage_property'	=> -1
};


sub new {
    my ($class, %args) = @_;
    my $self = $class->SUPER::new(%args );
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
sub get_study_stage_properties {
    my $self = shift;
    return IPD::IPDObject::StudyStageProperty->get_study_stage_properties_by_study_stage( $self->{'id'} );
}

#Class and instance methods using this class

=item IPD::IPDObject::StudyStage->get_study_stage($id)

B<Description:> Retrieves study stage information from a single study stage ID from IPD and stores into an object
B<Parameters:> The ID of the study stage to be queried.
B<Returns:> An IPD::IPDObject::StudyStage object

=cut

sub get_study_stage {
    my $class = shift;
    die ("Must call 'get' methods in the form of ", ref($class), "->method() when returning an object of this type: $!\n") if ref $class;
    my $id = shift;

    $extension = "study_stages/$id.xml";
    my $elem = $class->SUPER::XML2Object(IPD::Net->send_request(IPD::Client->get_client, $extension, "GET"));
    my $ss = $class->new(%{$elem});
    return $ss;
}

=item IPD::IPDObject::StudyStage->get_study_stages_by_study($id)

B<Description:> Retrieves all study stage information associated with a single study ID from IPD and stores into an array of objects
B<Parameters:> The ID of the study to be queried.
B<Returns:> An array of IPD::IPDObject::StudyStage objects

=cut

sub get_study_stages_by_study {
    my $class = shift;
    die ("Must call 'get' methods in the form of ", ref($class), "->method() when returning an object of this type: $!\n") if ref $class;
    my $id = shift;
    my @selves;

    $extension = "study_stages.xml?study_id=$id";
    my $elem = $class->SUPER::XML2Object(IPD::Net->send_request(IPD::Client->get_client, $extension, "GET"));
    if (ref($elem->{'study_stage'}) eq 'HASH') {
	push @selves, IPD::IPDObject::StudyStage->new(%{$elem->{'study_stage'}});
    } else {
    	foreach my $ss ( @{$elem->{'study_stage'}} ) {
	    push @selves, $class->new(%{$ss}); 
        }
    }
    return \@selves;
}

=item $obj->create_study_stage()

B<Description:> Creates a new Study Stage entry in IPD based on the provided StudyStage object
B<Parameters:> An IPD::IPDObject::StudyStage object
B<Returns:> The same IPD::IPDObject::StudyStage object

=cut

sub create_study_stage {
    my $self = shift;
    die ("Must call 'create|save|delete' methods in the form of ", '$object->method()', " : $!\n") unless ref $self;

    foreach my $name ('study_id', 'type', 'status') {
	if (!exists(${$self}{$name})) {
	    die "ERROR:  study_id, type, and status elements are required. $!\n";	
	}
    }
    $extension = "study_stages.xml";
    my $xml =  ref($self)->SUPER::construct_xml($self, "study_stage");
    my $temp = ref($self)->SUPER::XML2Object(IPD::Net->send_request(IPD::Client->get_client, $extension, "POST", $xml));
    my $id = $temp->{'id'};
    $extension = "study_stages/$id.xml";
    my $elem = ref($self)->SUPER::XML2Object(IPD::Net->send_request(IPD::Client->get_client, $extension, "GET"));
    $self = ref($self)->new(%{$elem});
    return $self;
}

=item $obj->save_study_stage()

B<Description:> Updates a Study Stage entry in IPD based on the provided StudyStage object.  The object must have an 'id' attribute associated with an existing study stage in IPD
B<Parameters:> An IPD::IPDObject::StudyStage object
B<Returns:> The same IPD::IPDObject::StudyStage object

=cut

sub save_study_stage {
    my $self = shift;
    die ("Must call 'create|save|delete' methods in the form of ", '$object->method()', " : $!\n") unless ref $self;
    my $id;

    exists $self->{'id'} > 0 ? $id = $self->{'id'} : die "No Study Stage ID present.  Perhaps you wanted to use 'create_study_stage' method? $!\n";
    $extension = "study_stages/$id.xml";
    my $xml =  ref($self)->SUPER::construct_xml($self, "study_stage");
    IPD::Net->send_request(IPD::Client->get_client, $extension, "PUT", $xml);
    my $elem = ref($self)->SUPER::XML2Object(IPD::Net->send_request(IPD::Client->get_client, $extension, "GET"));
    $self = ref($self)->new(%{$elem});
    return $self;
}

=item $obj->delete_study_stage()

B<Description:> Deletes a Study Stage entry from IPD, based on the ID attribute found in the StudyStage object.  Also clears all attributes from the object.
B<Parameters:> An IPD::IPDObject::StudyStage object
B<Returns:> Nothing

=cut

sub delete_study_stage {
    my $self = shift;
    die ("Must call 'create|save|delete' methods in the form of ", '$object->method()', " : $!\n") unless ref $self;

    die "No Study Stage ID attribute present for delete command: $!\n" if (! exists $self->{'id'});
    my $id = $self->{'id'};
    $extension = "study_stages/$id.xml";
    IPD::Net->send_request(IPD::Client->get_client, $extension, "DELETE");
    $self = ();
    return;
}

1;
