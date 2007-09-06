package OBO::OBOTerm;

=head1 NAME

OBO::OBOTerm.pm

=head1 VERSION

1.0

=head1 SYNOPSIS

use OBO::OBOTerm;
my $obj = new OBO::OBOTerm($id, $name, $def, $is_obsolete);
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
use Data::Dumper;

my $logger = OBO::Logger::get_logger("Logger::OBO");



## This will ensure that the term tags are written to file in the correct order
my @TAG_ORDER = qw(id is_anonymous name namespace alt_id def comment subset synonym xref is_a intersection_of union_of disjoint_from relationship is_obsolete replaced_by consider);

my $HAS_MULTIPLE_VALUES = { 'alt_id'=> 1,
			    'xref'=> 1,
			    'is_a'=> 1,
			    'subset'=>1,
			    'synonym'=>1,
			    'intersection_of'=>1,
			    'union_of'=>1,
			    'disjoint_from'=>1,
			    'relationship'=>1,
			    'consider'=>1
			};

## These are the required term tags
my $REQUIRED_TAGS = { 'id' => 1,
		      'name' => 1 };

## To support nextXref()
my $xrefIndex=0;

## To support nextAltId()
my $altIdIndex=0;

## To support nextSynonym()
my $synonymIndex=0;

## To support nextExactSynonym()
my $exactSynonymIndex=0;

## To support nextBroadSynonym()
my $broadSynonymIndex=0;

## To support nextNarrowSynonym()
my $narrowSynonymIndex=0;

## To support nextRelationship()
my $relationshipIndex=0;

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
    
    if ((defined($is_obsolete)) && ($is_obsolete != 0)){
	$self->{'_is_obsolete'} = $is_obsolete;
    }
   
    bless $self, $class;    

    $self->_init(@_);

    return $self;
}

=item $self->_init()

B<Description:> Not being used at this time.

B<Parameters:> None at this time.

B<Returns:> Nothing at this time.

=cut 

sub _init {

    my $self = shift;

    ## xref
    $self->{'_xref_sorted'} = 0;

    $self->{'_sorted_xrefs'} = [];

    $self->{'_xref_ctr'} = 0;


    ## alt_id
    $self->{'_alt_ids_sorted'} = 0;

    $self->{'_sorted_alt_ids'} = [];

    $self->{'_alt_id_ctr'} = 0;

    ## synonym
    $self->{'_synonyms_sorted'} = 0;

    $self->{'_sorted_synonyms'} = [];

    $self->{'_synonym_ctr'} = 0;

    ## exact synonym
    $self->{'_exact_synonyms_sorted'} = 0;

    $self->{'_sorted_exact_synonyms'} = [];

    $self->{'_exact_synonym_ctr'} = 0;

    ## narrow synonym
    $self->{'_narrow_synonyms_sorted'} = 0;

    $self->{'_sorted_narrow_synonyms'} = [];

    $self->{'_narrow_synonym_ctr'} = 0;

    ## broad synonym
    $self->{'_broad_synonyms_sorted'} = 0;

    $self->{'_sorted_broad_synonyms'} = [];

    $self->{'_broad_synonym_ctr'} = 0;

    ## relationship
    $self->{'_relationship_ctr'} = 0;

}

=item DESTROY

B<Description:> OBORecord class destructor

B<Parameters:> None

B<Returns:> None

=cut 

sub DESTROY  {
    my $self = shift;

}


=item $obj->getId()

B<Description:> Retrieves the ID attribute for the OBOTerm

B<Parameters:> None

B<Returns:> Scalar/string that is the ID assigned to the OBOTerm

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
	$logger->logdie("OBOTerm with no id!");
    }
}

=item $obj->setId()

B<Description:> Set the id for the OBOTerm

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

B<Description:> Retrieves the name for the OBOTerm

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
	$logger->logdie("OBOTerm with no name!");
    }
}

=item $obj->setName()

B<Description:> Set the name for the OBOTerm

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

B<Description:> Retrieves the namespace for the OBOTerm

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

=item $obj->hasNamespace()

B<Description:> Checks whether the OBOTerm has a namespace value

B<Parameters:> None

B<Returns:> 

  0 - false (scalar)
  1 - true  (scalar)

=cut 

sub hasNamespace {

    my ($self) = shift;

    if (!defined($self)){
	$logger->logdie("self was not defined");
    }
    
    if ( exists $self->{'_namespace'} ) {
	return 1;
    }

    return 0;
}

=item $obj->setNamespace()

B<Description:> Set the namespace for the OBOTerm

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

B<Description:> Retrieves the def for the OBOTerm

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
	$logger->warn("def is not defined for OBOTerm with id '$self->{'_id'}'");
    }
}

=item $obj->setDef()

B<Description:> Set the def for the OBOTerm

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

B<Description:> Retrieves the comment for the OBOTerm

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
	$logger->warn("comment is not defined for OBOTerm with id '$self->{'_id'}'");
    }
}

=item $obj->setComment()

B<Description:> Set the comment for the OBOTerm

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

B<Description:> Add alt_id for the OBOTerm

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

    $self->{'_alt_id_ctr'}++;
}

=item $obj->hasAltId()

B<Description:> Check whether this OBOTerm has any alt_id attributes

B<Parameters:> none

B<Returns:> 
 
  0 - false, scalar
  1 - true, scalar

=cut 

sub hasAltId {

    my ($self) = shift;

    if (!defined($self)){
	$logger->logdie("self was not defined");
    }

    if ((exists $self->{'_alt_id_ctr'}) && ($self->{'_alt_id_ctr'} > 0)){
	return 1;
    }
    return 0;
}

=item $obj->addExactSynonym()

B<Description:> Add exact synonym for the OBOTerm

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

=item $obj->addSynonym()

B<Description:> Add synonym for the OBOTerm

B<Parameters:> $synonym - scalar

B<Returns:> none

=cut 

sub addSynonym {

    my ($self) = shift;
    my ($synonym) = @_;

    if (!defined($self)){
	$logger->logdie("self was not defined");
    }

    if (!defined($synonym)){
	$logger->logdie("synonym not defined");
    }
    
    $self->{'_synonym'}->{$synonym}++;
}

=item $obj->hasSynonym()

B<Description:> Check whether this OBOTerm has any synonym attributes

B<Parameters:> none

B<Returns:> 

  0 - false, scalar
  1 - true, scalar

=cut 

sub hasSynonym {

    my ($self) = shift;

    if (!defined($self)){
	$logger->logdie("self was not defined");
    }

    if ((exists $self->{'_synonym_ctr'}) && ($self->{'_synonym_ctr'} > 0)){
	return 1;
    }
    return 0;
}

=item $obj->addNarrowSynonym()

B<Description:> Add narrow synonym for the OBOTerm

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

=item $obj->hasNarrowSynonym()

B<Description:> Check whether this OBOTerm has any NARRROW synonym attributes

B<Parameters:> none

B<Returns:> 

  0 - false, scalar
  1 - true, scalar

=cut 

sub hasNarrowSynonym {

    my ($self) = shift;

    if (!defined($self)){
	$logger->logdie("self was not defined");
    }

    if ((exists $self->{'_narrow_synonym_ctr'}) && ($self->{'_narrow_synonym_ctr'} > 0)){
	return 1;
    }
    return 0;
}

=item $obj->addBroadSynonym()

B<Description:> Add broad synonym for the OBOTerm

B<Parameters:> $synonym - scalar

B<Returns:> none

=cut 

sub addBroadSynonym {

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

=item $obj->hasBroardSynonym()

B<Description:> Check whether this OBOTerm has any BROAD synonym attributes

B<Parameters:> none

B<Returns:> 

  0 - false, scalar
  1 - true, scalar

=cut 

sub hasBroadSynonym {

    my ($self) = shift;

    if (!defined($self)){
	$logger->logdie("self was not defined");
    }

    if ((exists $self->{'_broad_synonym_ctr'}) && ($self->{'_broad_synonym_ctr'} > 0)){
	return 1;
    }
    return 0;
}

=item $obj->addXref()

B<Description:> Add xref for the OBOTerm

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

    $self->{'_xref_ctr'}++;
}

=item $obj->hasXref()

B<Description:> Check whether this OBOTerm has any xef attributes

B<Parameters:> none

B<Returns:> 

  0 - false, scalar
  1 - true, scalar

=cut 

sub hasXref {

    my ($self) = shift;

    if (!defined($self)){
	$logger->logdie("self was not defined");
    }

    if ((exists $self->{'_xref_ctr'}) && ($self->{'_xref_ctr'} > 0)){
	return 1;
    }
    return 0;
}


=item $obj->addIsA($id)

B<Description:> Add is_a to the OBOTerm object

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

B<Description:> Add a relationship to the OBOTerm object

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

    if ($reltype eq 'is_a'){
	$self->addIsA($value);
    }
    else {
	$self->{'_relationship'}->{$reltype}->{$value}++;
    }

    $self->{'_relationship_ctr'}++;
}

=item $obj->hasRelationship()

B<Description:> Checks whether the OBOTerm has at least one relationship value

B<Parameters:> None

B<Returns:> 

  0 - false (scalar)
  1 - true  (scalar)

=cut 

sub hasRelationship {

    my ($self) = shift;

    if (!defined($self)){
	$logger->logdie("self was not defined");
    }
    
    if (( exists $self->{'_relationship_ctr'} ) && ($self->{'_relationship_ctr'} > 0)){ 
	return 1;
    }

    return 0;
}

=item $obj->nextRelationship()

B<Description:> Return the next relationship tuple

B<Parameters:> None

B<Returns:> 

  Reference to two-element array
  element 1 == $id - id of related Term (scalar)
  element 2 == $reltype - relationship type (scalar)

=cut 

sub nextRelationship {

    my ($self) = shift;

    if (!defined($self)){
	$logger->logdie("self was not defined");
    }


    if ($self->{'_relationships_sorted'} == 0 ) {

	## If the relationhships aren't already sorted, do so now.
	foreach my $reltype (sort keys %{$self->{'_relationship'}} ) {
	    foreach my $value (sort keys %{$self->{'_relationship'}->{$reltype}} ) {
		push ( @{$self->{'_sorted_relationships'}}, [$value, $reltype] );	
	    }
	}

	## The list has been sorted once, let's not do it again.
	$self->{'_relationships_sorted'} = 1;
    }

    if ( $relationshipIndex < $self->{'_relationship_ctr'} ){
	return $self->{'_sorted_relationships'}->[$relationshipIndex++];
    }
    
    return undef;
}

=item $obj->setComment($comment)

B<Description:> Set the comment for the OBOTerm object

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

=item $obj->hasComment()

B<Description:> Checks whether the Term has a comment

B<Parameters:> None

B<Returns:> 

 0 - false (scalar)
 1 - true  (scalar)

=cut 

sub hasComment {

    my ($self) = shift;

    if (!defined($self)){
	$logger->logdie("self was not defined");
    }

    if ((exists $self->{'_comment'}) && (defined($self->{'_comment'}))){
	return 1;
    }
    return 0;
}


=item $obj->getIsObsolete()

B<Description:> Returns the is_obsolete value for the OBOTerm

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

=item $obj->hasIsObsolete()

B<Description:> Determine whether a value exists for the is_obsolete attribute

B<Parameters:> None

B<Returns:>  

 0 - false
 1 - true

=cut 

sub hasIsObsolete {

    my ($self) = shift;

    if ( exists $self->{'_is_obsolete'} ) {
	return  1;
    }
    return 0;

}

=item $obj->setIsObsolete()

B<Description:> Sets the is_obsolete value for the OBOTerm

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

    print $fh "[Term]\n";

    foreach my $tag (@TAG_ORDER){
	
	my $privatetag = '_' . $tag;
	
	if ($tag eq 'synonym'){
	    $self->writeSynonyms($fh);
	}

	if (exists $self->{$privatetag}){
	    
	    if (exists $HAS_MULTIPLE_VALUES->{$tag}){
		
		if ($tag eq 'is_a'){
		    $self->writeIsA($fh);
		}
		elsif ($tag eq 'relationship'){
		    $self->writeRelationships($fh);
		}
		else {
		    
		    foreach my $val (sort keys %{$self->{$privatetag}}){
			print $fh "$tag: $val\n";
		    }
		}
	    }
	    else {
		if ($tag eq 'is_obsolete'){
		    if ((exists $self->{$privatetag}) && ($self->{$privatetag} != 0)){
			print $fh "$tag: true\n";
		    }
		}
		else {
		    print $fh "$tag: $self->{$privatetag}\n";
		}
	    }
	}
	else {
	    if ( exists $REQUIRED_TAGS->{$tag}){
		$logger->logdie("required Term tag '$tag' was not defined for OBOTerm object with id '$self->{'_id'}'");
	    }
	}
    }
}

=item $obj->writeSynonyms($fh)

B<Description:> Method that writes the synonyms for the OBO record to the filehandle fh.

B<Parameters:>

 $fh     - filehandle for output file

B<Returns:>  None

=cut 

sub writeSynonyms {

    my ($self) = shift;
    my ($fh) = @_;

    if (!defined($self)){
	$logger->logdie("self was not defined");
    }
    if (!defined($fh)){
	$logger->logdie("fh was not defined");
    }

    my $synonymLookup = { '_synonym' => 'SYNONYM',
			  '_exact_synonym' => 'EXACT',
			  '_narrow_synonym' => 'NARROW',
			  '_broad_synonym' => 'BROAD',
			  '_related_synonym' => 'RELATED' };

    foreach my $tag ( keys %{$synonymLookup} ){
	
	if (exists $self->{$tag}){
	    
	    foreach my $value (sort keys %{$self->{$tag}}){
		
		print $fh "synonym: \"$value\" $synonymLookup->{$tag} []\n";
	    }
	}
    }
}

=item $obj->writeIsA($fh)

B<Description:> Method that writes the is_a relationships for the OBO record to the filehandle fh.

B<Parameters:>

 $fh     - filehandle for output file

B<Returns:>  None

=cut 

sub writeIsA {

    my ($self) = shift;
    my ($fh) = @_;

    if (!defined($self)){
	$logger->logdie("self was not defined");
    }
    if (!defined($fh)){
	$logger->logdie("fh was not defined");
    }

    foreach my $isa ( sort keys %{$self->{'_is_a'}} ){
	
	my $name = $self->{'_builder'}->getRecordNameById($isa);
	if (defined($name)){
		
	    print $fh "is_a: $isa ! $name\n";
	}
	else {
	    print $fh "is_a: $isa\n";
	}
    }
}

=item $obj->writeRelationships($fh)

B<Description:> Method that writes the relationships for the OBO record to the filehandle fh.

B<Parameters:>

 $fh     - filehandle for output file

B<Returns:>  None

=cut 

sub writeRelationships {

    my ($self) = shift;
    my ($fh) = @_;

    if (!defined($self)){
	$logger->logdie("self was not defined");
    }
    if (!defined($fh)){
	$logger->logdie("fh was not defined");
    }

    foreach my $reltype ( sort keys %{$self->{'_relationship'}} ){
	foreach my $value ( sort keys %{$self->{'_relationship'}->{$reltype}} ){
	
	    my $name = $self->{'_builder'}->getRecordNameById($value);
	    if (defined($name)){
		
		print $fh "relationship: $reltype $value ! $name\n";
	    }
	    else {
		print $fh "relationship: $reltype $value\n";
	    }
	}
    }
}

=item $obj->nextXref()

B<Description:> Iteratively returns each xref

B<Parameters:> None

B<Returns:> $xref - scalar

=cut

sub nextXref {

    my ($self) = shift;

    if ((exists $self->{'_xref_ctr'}) && ($self->{'_xref_ctr'} > 0 )){
	
	if ($self->{'_xref_sorted'} == 0 ) {
	    
	    ## If the xrefs aren't already sorted, do so now.
	    foreach my $xref (sort keys %{$self->{'_xref'}} ) {
		push ( @{$self->{'_sorted_xrefs'}}, $xref );	
	    }
	    
	    ## The list has been sorted once, let's not do it again.
	    $self->{'_xref_sorted'} = 1;
	}
	
	if ( $xrefIndex < $self->{'_xref_ctr'} ){
	    my $xref = $self->{'_sorted_xrefs'}->[$xrefIndex++];
#	    print Dumper $xref;
	    if (defined($xref)){
		return $xref;
	    }
	    else {
		$logger->logdie("xref was not defined for xrefIndex '$xrefIndex'");
	    }
	}
    }
    else {
	$logger->logdie("Attempting to retrieve the next xref for id '$self->{'_id'}'");
    }
    return undef;
}

=item $obj->nextAltId()

B<Description:> Iteratively returns each alt_id

B<Parameters:> None

B<Returns:> $alt_id - scalar

=cut

sub nextAltId {

    my ($self) = shift;

    if ($self->{'_alt_ids_sorted'} == 0 ) {

	## If the alt_ids aren't already sorted, do so now.
	foreach my $alt_id (sort keys %{$self->{'_alt_id'}} ) {
	    push ( @{$self->{'_sorted_alt_ids'}}, $alt_id );	
	}

	## The list has been sorted once, let's not do it again.
	$self->{'_alt_ids_sorted'} = 1;
    }

    if ( $altIdIndex < $self->{'_alt_id_ctr'} ){
	my $alt_id = $self->{'_sorted_alt_ids'}->[$altIdIndex++];
	if (defined($alt_id)){
	    return $alt_id;
	}
	else {
	    $logger->logdie("alt_id was not defined for altIdIndex '$altIdIndex'");
	}
    }
    
    return undef;
}

=item $obj->nextSynonym()

B<Description:> Iteratively returns each synonym

B<Parameters:> None

B<Returns:> $synonym - scalar

=cut

sub nextSynonym {

    my ($self) = shift;

    if ($self->{'_synonyms_sorted'} == 0 ) {

	## If the synonyms aren't already sorted, do so now.
	foreach my $synonym (sort keys %{$self->{'_synonym'}} ) {
	    push ( @{$self->{'_sorted_synonyms'}}, $synonym );	
	}

	## The list has been sorted once, let's not do it again.
	$self->{'_synonyms_sorted'} = 1;
    }

    if ( $synonymIndex < $self->{'_synonym_ctr'} ){
	my $synonym = $self->{'_sorted_synonyms'}->[$synonymIndex++];
	if (defined($synonym)){
	    return $synonym;
	}
	else {
	    $logger->logdie("synonym was not defined for synonymIndex '$synonymIndex'");
	}
    }
    
    return undef;
}

=item $obj->hasSynonym()

B<Description:> Determine whether there are any synonym values

B<Parameters:> None

B<Returns:> 

 0 - false
 1 - true

=cut

sub hasSynonym {

    my ($self) = shift;

    if (exists $self->{'_synonym_ctr'}){
	if ($self->{'_synonym_ctr'} > 0 ){
	    return 1;
	}
    }

    return 0;
}


=item $obj->nextExactSynonym()

B<Description:> Iteratively returns each exact synonym

B<Parameters:> None

B<Returns:> $exactSynonym - scalar

=cut

sub nextExactSynonym {

    my ($self) = shift;

    if ($self->{'_exact_synonyms_sorted'} == 0 ) {

	## If the exact synonyms aren't already sorted, do so now.
	foreach my $exactSynonym (sort keys %{$self->{'_exact_synonym'}} ) {
	    push ( @{$self->{'_sorted_exact_synonyms'}}, $exactSynonym );	
	}

	## The list has been sorted once, let's not do it again.
	$self->{'_exact_synonyms_sorted'} = 1;
    }

    if ( $exactSynonymIndex < $self->{'_exact_synonym_ctr'} ){
	my $exactSynonym = $self->{'_sorted_exact_synonyms'}->[$exactSynonymIndex++];
	if (defined($exactSynonym)){
	    return $exactSynonym;
	}
	else {
	    $logger->logdie("exact synonym was not defined for exactSynonymIndex '$exactSynonymIndex'");
	}
    }
    
    return undef;
}

=item $obj->hasExactSynonym()

B<Description:> Determine whether there are any exact_synonym values

B<Parameters:> None

B<Returns:> 

 0 - false
 1 - true

=cut

sub hasExactSynonym {

    my ($self) = shift;

    if (exists $self->{'_exact_synonym_ctr'}){
	if ($self->{'_exact_synonym_ctr'} > 0 ){
	    return 1;
	}
    }

    return 0;
}


=item $obj->nextNarrowSynonym()

B<Description:> Iteratively returns each narrow synonym

B<Parameters:> None

B<Returns:> $narrowSynonym - scalar

=cut

sub nextNarrowSynonym {

    my ($self) = shift;

    if ($self->{'_narrow_synonyms_sorted'} == 0 ) {

	## If the narrow synonyms aren't already sorted, do so now.
	foreach my $narrowSynonym (sort keys %{$self->{'_narrow_synonym'}} ) {
	    push ( @{$self->{'_sorted_narrow_synonyms'}}, $narrowSynonym );	
	}

	## The list has been sorted once, let's not do it again.
	$self->{'_narrow_synonyms_sorted'} = 1;
    }

    if ( $narrowSynonymIndex < $self->{'_narrow_synonym_ctr'} ){
	my $narrowSynonym = $self->{'_sorted_narrow_synonyms'}->[$narrowSynonymIndex++];
	if (defined($narrowSynonym)){
	    return $narrowSynonym;
	}
	else {
	    $logger->logdie("narrow synonym was not defined for narrowSynonymIndex '$narrowSynonymIndex'");
	}
    }
    
    return undef;
}

=item $obj->hasNarrowSynonym()

B<Description:> Determine whether there are any narrow_synonym values

B<Parameters:> None

B<Returns:> 

 0 - false
 1 - true

=cut

sub hasNarrowSynonym {

    my ($self) = shift;

    if (exists $self->{'_narrow_synonym_ctr'}){
	if ($self->{'_narrow_synonym_ctr'} > 0 ){
	    return 1;
	}
    }

    return 0;
}


=item $obj->nextBroadSynonym()

B<Description:> Iteratively returns each broad synonym

B<Parameters:> None

B<Returns:> $broadSynonym - scalar

=cut

sub nextBroadSynonym {

    my ($self) = shift;

    if ($self->{'_broad_synonyms_sorted'} == 0 ) {

	## If the broad synonyms aren't already sorted, do so now.
	foreach my $broadSynonym (sort keys %{$self->{'_broad_synonym'}} ) {
	    push ( @{$self->{'_sorted_broad_synonyms'}}, $broadSynonym );	
	}

	## The list has been sorted once, let's not do it again.
	$self->{'_broad_synonyms_sorted'} = 1;
    }

    if ( $broadSynonymIndex < $self->{'_broad_synonym_ctr'} ){
	my $broadSynonym = $self->{'_sorted_broad_synonyms'}->[$broadSynonymIndex++];
	if (defined($broadSynonym)){
	    return $broadSynonym;
	}
	else {
	    $logger->logdie("broad synonym was not defined for broadSynonymIndex '$broadSynonymIndex'");
	}
    }
    
    return undef;
}

=item $obj->hasBroadSynonym()

B<Description:> Determine whether there are any broad_synonym values

B<Parameters:> None

B<Returns:> 

 0 - false
 1 - true

=cut

sub hasBroadSynonym {

    my ($self) = shift;

    if (exists $self->{'_broad_synonym_ctr'}){
	if ($self->{'_broad_synonym_ctr'} > 0 ){
	    return 1;
	}
    }

    return 0;
}

1;
