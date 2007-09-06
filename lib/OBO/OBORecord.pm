package OBO::OBORecord;

=head1 NAME

OBO::OBORecord.pm

=head1 VERSION

1.0

=head1 SYNOPSIS

use OBORecord;
my $obj = new OBO::OBORecord($id, $name, $def, $is_obsolete);
$obj->writeRecord($fileHandle1);

=head1 AUTHOR

Jay Sundaram
sundaram@jcvi.org

=head1 METHODS


=over 4

=cut

use strict;
use OBO::OBOBuilder;
use OBO::Logger;

my $logger = OBO::Logger::get_logger("Logger::OBO");



## This will ensure that the term tags are written to file in the correct order
my @TERM_TAG_ORDER = qw(id is_anonymous name namespace alt_id def comment subset synonym xref is_a intersection_of union_of disjoint_from relationship is_obsolete replaced_by consider);

my $HAS_MULTIPLE_VALUES = { 'alt_id'=> 1,
			    'xref'=> 1,
			    'is_a'=> 1,
			    'subset'=>1,
			    'synonym'=>1,
			    'is_a'=>1,
			    'intersection_of'=>1,
			    'union_of'=>1,
			    'disjoint_from'=>1,
			    'relationship'=>1,
			    'consider'=>1
			};

## This will ensure that the typedef tags are written to file in the correct order
my @TYPEDEF_TAG_ORDER = qw(id is_anonymous name namespace alt_id def comment subset synonym xref domain range is_anti_symmetric is_cyclic is_reflexive is_symmetric is_transitive is_a inverse_of transitive_over relationship is_obsolete replaced_by consider);

## This will ensure that the instance tags are written to file in the correct order
my @INSTANCE_TAG_ORDER = qw(id is_anonymous name namespace alt_id comment synonym xref instance_of property_value is_obsolete replaced_by consider);

## These are the required term tags
my $REQUIRED_TERM_TAGS = { 'id' => 1,
			   'name' => 1 };

## These are the required instance tags
my $REQUIRED_INSTANCE_TAGS = { 'id' => 1,
			       'name' => 1,
			       'instance_of' => 1};




## Class variable for keeping track of the number of OBO records
my $oboRecordCtr = 0;

=item new($id, $name, $def, $is_obsolete)

B<Description:> Instantiate OBORecord object

B<Parameters:> 

 $id          - id
 $name        - name
 $def         - def
 $is_obsolete - is_obsolete

B<Returns:> reference to the OBORecord

=cut

sub new  {

    my $class = shift;

    my ($id, $name, $def, $is_obsolete) = @_;

    my $self = {};

    if (!defined($id)){
	$logger->logdie("id was not defined");
    }
    else {
	$self->{'_id'} = $id;
    }

    if (defined($name)){
	$self->{'_name'} = $name;
    }
    else {
	$logger->logdie("name was not defined");
    }

    if (defined($def)){
	$self->{'_def'} = $def;
    }
    
    if (defined($is_obsolete)){
	$self->{'_is_obsolete'} = $is_obsolete;
    }
  

    bless $self, $class;    

    return $self;
}

=item $self->_init()

B<Description:> Not being used at this time.

B<Parameters:> None at this time.

B<Returns:> Nothing at this time.

=cut 

sub _init {

    my $self = shift;
    
    ## Column 9: attributes
    $self->{'_attrs'} = {};

    $self->{'_obo_record_type'} = 'term';
}

=item DESTROY

B<Description:> OBORecord class destructor

B<Parameters:> None

B<Returns:> None

=cut 

sub DESTROY  {
    my $self = shift;

}

=item $obj->addAttribute($key, $value)

B<Description:> Add named attributes to the GFFRecord object

B<Parameters:>

 $key   - name of the attribute
 $value - value to add to the record

B<Returns:>  None

=cut 

sub addAttribute {

    my ($self) = shift;
    my ($key, $value) = @_;

    if (!defined($key)){
	$logger->logdie("key was not defined");
    }
    if (!defined($value)){
	$logger->logdie("value was not defined");
    }

    ## URL encoding for compliance with OBO spec
    $value =~ s/\\n/\n/g;
    $value =~ s/([\f\n\r\t;=%&,[:cntrl:]])/sprintf("%%%02X", ord($1))/seg;

    push(@{$self->{'_attrs'}->{$key}}, $value);
}


=item $obj->getParent()

B<Description:> Retrieves the Parent attribute for the GFFRecord object

B<Parameters:> None

B<Returns:>  Scalar/string that is the value assigned to the Parent attribute for the OBORecord

=cut 

sub getParent {

    my ($self) = shift;

    if (!defined($self)){
	$logger->logdie("self was not defined");
    }
    if ( exists $self->{'_parent'} ) {
	return $self->{'_parent'};
    }
    else {
	if ($logger->is_debug()){
	    $logger->debug("_parent does not exist for OBORecord with id '$self->{'_id'}'");
	}
    }
    
    return undef;
}

=item $obj->getId()

B<Description:> Retrieves the ID attribute for the OBORecord

B<Parameters:> None

B<Returns:> Scalar/string that is the ID assigned to the OBORecord

=cut 

sub getId {

    my ($self) = shift;

    if (!defined($self)){
	$logger->logdie("self was not defined");
    }
    
    if ( exists $self->{'_id'} ) {
	return  $self->{'_id'};
    }
    else {
	$logger->logdie("OBORecord with no id!");
    }
}

=item $obj->setId()

B<Description:> Set the id for the OBORecord

B<Parameters:> $id - scalar

B<Returns:> none

=cut 

sub setId {

    my ($self) = shift;
    my ($id) = @_;

    if (!defined($self)){
	$logger->logdie("self was not defined");
    }

    if (!defined($id)){
	$logger->logdie("id was not defined");
    }
    
    $self->{'_id'} = $id;
}

=item $obj->getName()

B<Description:> Retrieves the name for the OBORecord

B<Parameters:> None

B<Returns:> $name - scalar

=cut 

sub getName {

    my ($self) = shift;

    if (!defined($self)){
	$logger->logdie("self was not defined");
    }
    
    if ( exists $self->{'_name'} ) {
	return  $self->{'_name'};
    }
    else {
	$logger->logdie("OBORecord with no name!");
    }
}

=item $obj->setName()

B<Description:> Set the name for the OBORecord

B<Parameters:> $name - scalar

B<Returns:> none

=cut 

sub setName {

    my ($self) = shift;
    my ($name) = @_;

    if (!defined($self)){
	$logger->logdie("self was not defined");
    }

    if (!defined($name)){
	$logger->logdie("name was not defined");
    }
    
    $self->{'_name'} = $name;
}

=item $obj->getNamespace()

B<Description:> Retrieves the namespace for the OBORecord

B<Parameters:> None

B<Returns:> $namespace - scalar

=cut 

sub getNamespace {

    my ($self) = shift;

    if (!defined($self)){
	$logger->logdie("self was not defined");
    }
    
    if ( exists $self->{'_namespace'} ) {
	return  $self->{'_namespace'};
    }
    else {
	if (exists $self->{'_builder'}){
	    return $self->{'_builder'}->getDefaultNamespace();
	}
	else {
	    $logger->warn("namespace is not defined and could not retrieve default-namespace");
	}
    }
}

=item $obj->setNamespace()

B<Description:> Set the namespace for the OBORecord

B<Parameters:> $namespace - scalar

B<Returns:> none

=cut 

sub setNamespace {

    my ($self) = shift;
    my ($namespace) = @_;

    if (!defined($self)){
	$logger->logdie("self was not defined");
    }

    if (!defined($namespace)){
	$logger->logdie("namespace was not defined");
    }
    
    $self->{'_namespace'} = $namespace;
}

=item $obj->getDef()

B<Description:> Retrieves the def for the OBORecord

B<Parameters:> None

B<Returns:> $def - scalar

=cut 

sub getDef {

    my ($self) = shift;

    if (!defined($self)){
	$logger->logdie("self was not defined");
    }
    
    if ( exists $self->{'_def'} ) {
	return  $self->{'_def'};
    }
    else {
	$logger->warn("def is not defined for OBORecord with id '$self->{'_id'}'");
    }
}

=item $obj->setDef()

B<Description:> Set the def for the OBORecord

B<Parameters:> $def - scalar

B<Returns:> none

=cut 

sub setDef {

    my ($self) = shift;
    my ($def) = @_;

    if (!defined($self)){
	$logger->logdie("self was not defined");
    }

    if (!defined($def)){
	$logger->logdie("def was not defined");
    }
    
    $self->{'_def'} = $def;
}


=item $obj->getComment()

B<Description:> Retrieves the comment for the OBORecord

B<Parameters:> None

B<Returns:> $comment - scalar

=cut 

sub getComment {

    my ($self) = shift;

    if (!defined($self)){
	$logger->logdie("self was not defined");
    }
    
    if ( exists $self->{'_comment'} ) {
	return  $self->{'_comment'};
    }
    else {
	$logger->warn("comment is not defined for OBORecord with id '$self->{'_id'}'");
    }
}

=item $obj->setComment()

B<Description:> Set the comment for the OBORecord

B<Parameters:> $comment - scalar

B<Returns:> none

=cut 

sub setComment {

    my ($self) = shift;
    my ($comment) = @_;

    if (!defined($self)){
	$logger->logdie("self was not defined");
    }

    if (!defined($comment)){
	$logger->logdie("comment was not defined");
    }
    
    $self->{'_comment'} = $comment;
}

=item $obj->addAltId()

B<Description:> Add alt_id for the OBORecord

B<Parameters:> $alt_id - scalar

B<Returns:> none

=cut 

sub addAltId {

    my ($self) = shift;
    my ($alt_id) = @_;

    if (!defined($self)){
	$logger->logdie("self was not defined");
    }

    if (!defined($alt_id)){
	$logger->logdie("alt_id not defined");
    }
    
    $self->{'_alt_id'}->{$alt_id}++;
}

=item $obj->addExactSynonym()

B<Description:> Add exact synonym for the OBORecord

B<Parameters:> $synonym - scalar

B<Returns:> none

=cut 

sub addExactSynonym {

    my ($self) = shift;
    my ($synonym) = @_;

    if (!defined($self)){
	$logger->logdie("self was not defined");
    }

    if (!defined($synonym)){
	$logger->logdie("synonym not defined");
    }
    
    $self->{'_exact_synonym'}->{$synonym}++;
}

=item $obj->addRelatedSynonym()

B<Description:> Add related synonym for the OBORecord

B<Parameters:> $synonym - scalar

B<Returns:> none

=cut 

sub addRelatedSynonym {

    my ($self) = shift;
    my ($synonym) = @_;

    if (!defined($self)){
	$logger->logdie("self was not defined");
    }

    if (!defined($synonym)){
	$logger->logdie("synonym not defined");
    }
    
    $self->{'_related_synonym'}->{$synonym}++;
}

=item $obj->addNarrowSynonym()

B<Description:> Add narrow synonym for the OBORecord

B<Parameters:> $synonym - scalar

B<Returns:> none

=cut 

sub addNarrowSynonym {

    my ($self) = shift;
    my ($synonym) = @_;

    if (!defined($self)){
	$logger->logdie("self was not defined");
    }

    if (!defined($synonym)){
	$logger->logdie("synonym not defined");
    }
    
    $self->{'_narrow_synonym'}->{$synonym}++;
}

=item $obj->addBroadSynonym()

B<Description:> Add broad synonym for the OBORecord

B<Parameters:> $synonym - scalar

B<Returns:> none

=cut 

sub addBroardSynonym {

    my ($self) = shift;
    my ($synonym) = @_;

    if (!defined($self)){
	$logger->logdie("self was not defined");
    }

    if (!defined($synonym)){
	$logger->logdie("synonym not defined");
    }
    
    $self->{'_broad_synonym'}->{$synonym}++;
}


=item $obj->addXref()

B<Description:> Add xref for the OBORecord

B<Parameters:> $xref - scalar

B<Returns:> none

=cut 

sub addXref {

    my ($self) = shift;
    my ($xref) = @_;

    if (!defined($self)){
	$logger->logdie("self was not defined");
    }

    if (!defined($xref)){
	$logger->logdie("xref not defined");
    }
    
    $self->{'_xref'}->{$xref}++;
}

=item $obj->addIsA($id)

B<Description:> Add is_a to the OBORecord object

B<Parameters:> $id - scalar

B<Returns:> None

=cut 

sub addIsA {

    my ($self) = shift;
    my ($isa) = @_;

    if (!defined($self)){
	$logger->logdie("self was not defined");
    }
    if (!defined($isa)){
	$logger->logdie("isa was not defined");
    }

    $self->{'_is_a'}->{$isa}++;
}

=item $obj->addRelationship($reltype,$value)

B<Description:> Add a relationship to the OBORecord object

B<Parameters:> 

  $reltype - scalar
  $value   - scalar

B<Returns:> None

=cut 

sub addRelationship {

    my ($self) = shift;
    my ($reltype, $value) = @_;

    if (!defined($self)){
	$logger->logdie("self was not defined");
    }

    if (!defined($reltype)){
	$logger->logdie("reltype was not defined");
    }

    if (!defined($value)){
	$logger->logdie("value was not defined");
    }

    my $key = '_' . $reltype;

    $self->{$key}->{$value}++;
}

=item $obj->setComment($comment)

B<Description:> Set the comment for the OBORecord object

B<Parameters:> 

  $comment - scalar

B<Returns:> None

=cut 

sub setComment {

    my ($self) = shift;
    my ($comment) = @_;

    if (!defined($self)){
	$logger->logdie("self was not defined");
    }

    if (!defined($comment)){
	$logger->logdie("comment was not defined");
    }


    $self->{'_comment'} = $comment;
}


=item $obj->getIsObsolete()

B<Description:> Returns the is_obsolete value for the OBORecord

B<Parameters:> None

B<Returns:>  is_obsolete value (scalar)

=cut 

sub getIsObsolete {

    my ($self) = shift;

    if (!defined($self)){
	$logger->logdie("self was not defined");
    }
    
    if ( exists $self->{'_is_obsolete'} ) {
	return  $self->{'_is_obsolete'};
    }
    else {
	return 0;
    }
}

=item $obj->setIsObsolete()

B<Description:> Sets the is_obsolete value for the OBORecord

B<Parameters:> $is_obsolete (scalar)

B<Returns:>  none

=cut 

sub setIsObsolete {

    my $self = shift;
    my ($is_obsolete) = @_;

    if (!defined($self)){
	$logger->logdie("self was not defined");
    }

    if (!defined($is_obsolete)){
	$logger->logdie("is_obsolete was not defined");
    }
    else {
	$self->{'_is_obsolete'} = $is_obsolete;
    }
}

=item $obj->writeRecord($fh)

B<Description:> Method that writes the OBO record to the filehandle fh.

B<Parameters:>

 $fh     - filehandle for output file

B<Returns:>  None

=cut 

sub writeRecord {

    my ($self) = shift;
    my ($fh) = @_;

    if (!defined($self)){
	$logger->logdie("self was not defined");
    }
    if (!defined($fh)){
	$logger->logdie("fh was not defined");
    }


    if ($self->{'_obo_record_type'} eq 'term'){
	
	foreach my $tag (@TERM_TAG_ORDER){

	    my $privatetag = '_' . $tag;

	    if (exists $self->{$privatetag}){

		if (exists $HAS_MULTIPLE_VALUES->{$privatetag}){

		    foreach my $val (@{$self->{$privatetag}}){
			print $fh "$tag: $val\n";
		    }
		}
		else {
		    print $fh "$tag: $self->{$privatetag}\n";
		}
	    }
	    else {
		if ( exists $REQUIRED_TERM_TAGS->{$tag}){
		    $logger->logdie("required Term tag '$tag' was not defined for OBORecord with id '$self->{'_id'}'");
		}
	    }
	}
    }
    elsif ($self->{'_obo_record_type'} eq 'typedef'){

	foreach my $tag (@TYPEDEF_TAG_ORDER){

	    my $privatetag = '_' . $tag;

	    if (exists $self->{$privatetag}){

		if (exists $HAS_MULTIPLE_VALUES->{$privatetag}){

		    foreach my $val (@{$self->{$privatetag}}){
			print $fh "$tag: $val\n";
		    }
		}
		else {
		    print $fh "$tag: $self->{$privatetag}\n";
		}
	    }
	    else {
		if ( exists $REQUIRED_TERM_TAGS->{$tag}){
		    $logger->logdie("required Typedef tag '$tag' was not defined for OBORecord with id '$self->{'_id'}'");
		}
	    }
	}
    }
    else {
	$logger->logdie("This type of OBO record '$self->{'_obo_record_type'}' is not support");
    }
}

1;
