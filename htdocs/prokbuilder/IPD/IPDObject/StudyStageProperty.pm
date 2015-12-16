package IPD::IPDObject::StudyStageProperty;

=head1 NAME 

IPD::IPDObject::StudyStageProperty - Module for managing attributes of StudyStageProperty-type IPD objects

=cut

use IPD::Client;
use strict;
use warnings;
use base qw(IPD::IPDObject);

our $AUTOLOAD;

my $extension = "";

my $attributes = {
	'study_stage_id' => -1,
	'id' => -1,
	'type' => -1,
	'value' => -1,
	'updated_at' => -1,
	'created_at' => -1
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

#Class and instance methods using this class

=item IPD::IPDObject::StudyStageProperty->get_study_stage_property($id)

B<Description:> Retrieves study stage property information from a single study stage property ID from IPD and stores into an object
B<Parameters:> The ID of the study stage property to be queried.
B<Returns:> An IPD::IPDObject::StudyStageProperty object

=cut

sub get_study_stage_property {
    my $class = shift;
    die ("Must call 'get' methods in the form of ", ref($class), "->method() : $!\n") if ref $class;
    my $id = shift;

    $extension = "study_stage_properties/$id.xml";
    my $elem = $class->SUPER::XML2Object(IPD::Net->send_request(IPD::Client->get_client, $extension, "GET"));
    my $ssp = $class->new(%{$elem});
    return $ssp;
}

=item IPD::IPDObject::StudyStageProperty>get_study_stage_properties_by_study_stage($id)

B<Description:> Retrieves all study stage property information associated with a single study stage ID from IPD and stores into an array of objects
B<Parameters:> The ID of the study stage to be queried.
B<Returns:> An array of IPD::IPDObject::StudyStageProperty objects

=cut

sub get_study_stage_properties_by_study_stage {
    my $class = shift;
    die ("Must call 'get' methods in the form of ", ref($class), "->method() : $!\n") if ref $class;
    my $id = shift;
    my @selves;
    my @ssp;

    $extension = "study_stages/$id.xml";
    my $elem = $class->SUPER::XML2Object(IPD::Net->send_request(IPD::Client->get_client, $extension, "GET"));
    if (ref($elem->{'study_stage_properties'}->{'study_stage_property'}) eq 'HASH') {
	push @ssp, $elem->{'study_stage_properties'}->{'study_stage_property'}->{'id'};
    } else {
    	foreach my $prop ( @{$elem->{'study_stage_properties'}->{'study_stage_property'}} ) {
	    push @ssp, $prop->{'id'};
        }
    }
    foreach my $ssp_id (@ssp) {
	push @selves, $class->get_study_stage_property($ssp_id);
    }
    return \@selves;
}


=item $obj->create_study_stage_property()

B<Description:> Creates a new Study Stage Property entry in IPD based on the provided StudyStageProperty object
B<Parameters:> An IPD::IPDObject::StudyStageProperty object
B<Returns:> The same IPD::IPDObject::StudyStageProperty object

=cut

sub create_study_stage_property {
    my $self = shift;
    die ("Must call 'create|save|delete' methods in the form of ", '$object->method()', " : $!\n") unless ref $self;

    foreach my $name ('study_stage_id', 'type', 'value') {
	if (!exists(${$self}{$name})) {
	    die "ERROR:  study_stage_id, type, and value elements are required. $!\n";	
	}
    }
    $extension = "study_stage_properties.xml";
    my $xml =  ref($self)->SUPER::construct_xml($self, "study_stage_property");
    my $temp = ref($self)->SUPER::XML2Object(IPD::Net->send_request(IPD::Client->get_client, $extension, "POST", $xml));
    my $id = $temp->{'id'};
    $extension = "study_stage_properties/$id.xml";
    my $elem = ref($self)->SUPER::XML2Object(IPD::Net->send_request(IPD::Client->get_client, $extension, "GET"));
    $self = ref($self)->new(%{$elem});
    return $self;
}

=item $obj->save_study_stage_property()

B<Description:> Updates a Study Stage Property entry in IPD based on the provided StudyStageProperty object.  The object must have an 'id' attribute associated with an existing study stage property in IPD
B<Parameters:> An IPD::IPDObject::StudyStageProperty object
B<Returns:> The same IPD::IPDObject::StudyStageProperty object

=cut

sub save_study_stage_property {
    my $self = shift;
    die ("Must call 'create|save|delete' methods in the form of ", '$object->method()', " : $!\n") unless ref $self;
    my $id;

    exists $self->{'id'} ? $id = $self->{'id'} : die "No Study Stage Property ID present.  Perhaps you wanted to use 'create_study_stage_property' method? $!\n";
    $extension = "study_stage_properties/$id.xml";
    my $xml =  ref($self)->SUPER::construct_xml($self, "study_stage_property");
    IPD::Net->send_request(IPD::Client->get_client, $extension, "PUT", $xml);
    my $elem = ref($self)->SUPER::XML2Object(IPD::Net->send_request(IPD::Client->get_client, $extension, "GET"));
    $self = ref($self)->new(%{$elem});
    return $self;
}

=item $obj->delete_study_stage_property()

B<Description:> Deletes a Study Stage Property entry from IPD, based on the ID attribute found in the StudyStageProperty object.  Also clears all attributes from the object.
B<Parameters:> An IPD::IPDObject::StudyStageProperty object
B<Returns:> Nothing

=cut

sub delete_study_stage_property {
    my $self = shift;
    die ("Must call 'create|save|delete' methods in the form of ", '$object->method()', " : $!\n") unless ref $self;

    die "No Study Stage Property ID attribute present for delete command: $!\n" if (! exists $self->{'id'});
    my $id = $self->{'id'};
    $extension = "study_stage_properties/$id.xml";
    IPD::Net->send_request(IPD::Client->get_client, $extension, "DELETE");
    $self = ();
    return;
}


1;
