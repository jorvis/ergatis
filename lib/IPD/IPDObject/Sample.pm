package IPD::IPDObject::Sample;

=head1 NAME 

IPD::IPDObject::Sample - Module for managing attributes of Sample-type IPD objects

=cut

use IPD::Client;

use strict;
use warnings;
use base qw(IPD::IPDObject);

our $AUTOLOAD;

my $extension = "";

my $attributes = {
#common
	'name' => -1, 
	'sample_type' => -1,
	'id' => -1,
	'project_id' => -1,
	'study_id' => -1,
	'geographic_location' => -1,
	'type' => -1,
	'updated_at' => -1,
	'created_at' => -1,
	'lims_sample_ids' => -1,
	'habitat' => -1,
	'source_materials' => -1,
	'isolation_and_growth_condition' => -1,
	'nucleic_acid_preparation' => -1,
	'assembly_method' => -1,
	'sequencing_method' => -1,
	'sample_acquisition_date' => -1,
	'locus_tag_prefix' => -1,
	'lims_sample_ids' => -1,
#metagenomics
	'collection_time' => -1,
	'volume_of_sample' => -1,
	'sampling_strategy' => -1,
	'biomaterial_treatment' => -1,
	'metagenomic_sample_type' => -1,
#genomics
	'subspecific_genetic_lineage' => -1,
	'ploidy' => -1,
	'estimated_size_before_sequencing' => -1,
	'number_of_replicons' => -1,
	'extrachromosomal_elements' => -1,
	'pathogenicity' => -1,
	'biotic_relationship' => -1,
	'specific_host' => -1,
	'host_specificity' => -1,
	'health_or_disease_status' => -1,
	'trophic_level' => -1,
	'propagation' => -1,
	'encoded_trait' => -1,
	'relationship_to_oxygen' => -1,
	'organism_name' => -1,
	'organism_type' => -1,
	'strain' => -1,
	'common_name' => -1,
	'taxon_id' => -1
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
sub get_lims_sample_ids {
    my ($self) = @_;
    my @lims;
    if (ref($self->{'lims_sample_ids'}->{'lims_sample_id'}) eq 'HASH') {
	push @lims, $self->{'lims_sample_ids'}->{'lims_sample_id'};
    } else {
    	foreach my $sample ( @{$self->{'lims_sample_ids'}->{'lims_sample_id'}} ) {
	    push @lims, $sample;
        }
    }
    return \@lims;	#will return array-ref
}

#Methods we do not want to access with AUTOLOAD
sub set_id {
    my $self = shift;
    my $value = shift;
    die ("Not allowed to change the ID: $!\n") if (defined $self -> {'id'});
    $self -> {'id'} = $value;
}

sub set_sample_type {
    my $self = shift;
    my $value = shift;
    die ("Not allowed to change the sample type: $!\n") if (defined $self->{'sample_type'});
    $self->{'sample_type'} = $value;
}

#Class and instance methods using this class

=item IPD::IPDObject::Sample->get_sample($id)

B<Description:> Retrieves sample information from a single sample ID from IPD and stores into an object
B<Parameters:> The ID of the sample to be queried.
B<Returns:> An IPD::IPDObject::Sample object

=cut

sub get_sample {
    my $class = shift;
    die ("Must call 'get' methods in the form of ", ref($class), "->method() : $!\n") if ref $class;
    my $id = shift;

    $extension = "samples/$id.xml";
    my $elem = $class->SUPER::XML2Object(IPD::Net->send_request(IPD::Client->get_client, $extension, "GET"));
    my $sample = $class->new(%{$elem});
    return $sample;
}

=item IPD::IPDObject::Sample->get_all_samples

B<Description:> Retrieves all sample information from IPD and stores into an array of objects
B<Parameters:> None.
B<Returns:> An array of IPD::IPDObject::Sample objects
B<Note:>  Because of the sheer amount of samples currently present in IPD, this command may give a "time out" error due to taking too long (the time limit set is 999 seconds).  If this error occurs, contact me (sadkins@som.umaryland.edu) and I will make the limit longer

=cut

sub get_all_samples {
    my $class = shift;
    die ("Must call 'get' methods in the form of ", ref($class), "->method() : $!\n") if ref $class;
    my @selves; 

    $extension = "samples.xml";
    my $elem = $class->SUPER::XML2Object(IPD::Net->send_request(IPD::Client->get_client, $extension, "GET"));
   if (ref($elem->{'sample'}) eq 'HASH') {
	push @selves, $class->new(%{$elem->{'sample'}});
    } else {
    	foreach my $sample ( @{$elem->{'sample'}} ) {
	    push @selves, $class->new(%{$sample}); 
        }
    }
    return \@selves;
}

=item IPD::IPDObject::Sample->get_samples_by_project($id)

B<Description:> Retrieves all sample information associated with a project ID from IPD and stores into an array of objects
B<Parameters:> A project ID.
B<Returns:> An array of IPD::IPDObject::Sample objects

=cut

sub get_samples_by_project {
    my $class = shift;
    die ("Must call 'get' methods in the form of ", ref($class), "->method() : $!\n") if ref $class;
    my $id = shift;
    my @selves;

    $extension = "samples.xml?project_id=$id";
    my $elem = $class->SUPER::XML2Object(IPD::Net->send_request(IPD::Client->get_client, $extension, "GET"));
    if (ref($elem->{'sample'}) eq 'HASH') {
	push @selves, $class->new(%{$elem->{'sample'}});
    } else {
    	foreach my $sample ( @{$elem->{'sample'}} ) {
	    push @selves, $class->new(%{$sample}); 
        }
    }
    return \@selves;
}

=item IPD::IPDObject::Sample->get_sample_by_study($id)

B<Description:> Retrieves a sample associated with a study ID from IPD and stores into an array of objects
B<Parameters:> A study ID.
B<Returns:> An IPD::IPDObject::Sample object

=cut

sub get_sample_by_study {
    my $class = shift;
    die ("Must call 'get' methods in the form of ", ref($class), "->method() : $!\n") if ref $class;
    my $id = shift;
    my $sample;

    $extension = "studies/$id.xml";
    my $elem = $class->SUPER::XML2Object(IPD::Net->send_request(IPD::Client->get_client, $extension, "GET"));
    $sample = $class->get_sample($elem->{'sample'}->{'id'});
    return $sample;
}

1;
