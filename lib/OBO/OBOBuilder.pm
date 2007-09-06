package OBO::OBOBuilder;
=head1 NAME

OBO::OBOBuilder.pm 

A class to facilitate the creation of OBO files.

=head1 VERSION

1.0

=head1 SYNOPSIS

use OBOBuilder;
my $obj = new OBO::OBOBuilder( filename => '/tmp/myfile.obo');
$obj->addRecord($oboRecord);
$obj->writeFile();

=head1 AUTHOR

Jay Sundaram
sundaram@jcvi.org

=head1 METHODS


=over 4

=cut

use strict;
use OBO::OBOTerm;
use OBO::OBOTypedef;
use OBO::Logger;



## This will ensure that the headers are written to file in the correct order
my @HEADER_ORDER = qw( format-version data-version date saved-by auto-generated-by import subsetdef synonymtypedef default-namespace remark);

## These are the required headers
my $REQUIRED_HEADERS = {'format-version' => 1};


my $logger = OBO::Logger::get_logger("Logger::OBO");


## Class variable for keeping track of the number of OBOTerm records
my $oboTermCtr = 0;

## Class variable for keeping track of the number of OBOTypedef records
my $oboTypedefCtr = 0;

## Instance variable for supporting the nextTerm() method
my $recordIndex = 0;

## Instance variable for supporting the nextTypedef() method
my $typedefRecordIndex = 0;

## Instance variable for supporting the nextTerm() method
my $sorted = 0;

## Instance variable for supporting the nextRecord() method
my $typedefs_sorted = 0;

## To ensure that the name-namespace tuples are not repeated
my $nameNamespaceLookup = {};

=item new()

B<Description:> Instantiate OBOBuilder object

B<Parameters:> 

%args

B<Returns:> Returns a reference to OBOBuilder

=cut

sub new  {

    my $class = shift;
    my $self = {};
    bless $self, $class;    

    $self->_init(@_);

    return $self;
}

=item $self->_init(%args)

B<Description:> Typical Perl init() method

B<Parameters:> 

%args

B<Returns:> None

=cut

sub _init {

    my $self = shift;
    my (%args) = @_;

    foreach my $key (keys %args ){
	$self->{"_$key"} = $args{$key};
    }

    ## Instance data members

    ## This lookup is keyed on OBO header name with value being the OBO header value.
    $self->{'_headers'} = {};

    ## This lookup is keyed on the id with value being a reference to the OBOTerm.
    ## Reference for every OBOTerm created by the OBOBuilder will be stored in
    ## this lookup.
    $self->{'_terms'} = {};

    ## This lookup is keyed on the id with value being a reference to the OBOTypedef.
    ## Reference for every OBOTerm created by the OBOBuilder will be stored in
    ## this lookup.
    $self->{'_typedefs'} = {};

    $self->{'_headers'}->{'format-version'} = 1.2;

}

=item DESTROY

B<Description:> OBOBuilder class destructor

B<Parameters:> None

B<Returns:> None

=cut

sub DESTROY  {
    my $self = shift;

}

=item $obj->addHeader($key, $value)

B<Description:> Add OBO file header

B<Parameters:> 

$key   - scalar/string name of the header
$value - scalar/string value of the header attribute

B<Returns:> None

=cut

sub addHeader {

    my ($self) = shift;
    my ($key, $value) = @_;

    if (!defined($key)){
	$logger->logdie("key was not defined");
    }
    if (!defined($value)){
	$logger->logdie("value was not defined");
    }

    $self->{'_headers'}->{$key} = $value;
}

=item $obj->doesHeaderExist($name)

B<Description:> Check whether a header already exists with given name

B<Parameters:> $name - scalar/string name of the header

B<Returns:> 

0 - scalar false
1 - scalar true

=cut

sub doesHeaderExist {

    my ($self) = shift;
    my ($name) = @_;

    if (!defined($name)){
	$logger->logdie("header name was not defined");
    }

    if (exists $self->{'_headers'}->{$name}){
	return 1;
    }
    return 0;
}

=item $obj->setFilename($filename)

B<Description:> Set the output OBO filename

B<Parameters:> $filename - scalar/string name of the output file

B<Returns:> None

=cut

sub setFilename {

    my ($self) = shift;
    my ($filename) = @_;

    if (!defined($filename)){
	$logger->logdie("output filename was not defined");
    }

    $self->{'_filename'} = $filename;
}

=item $obj->getFilename($filename)

B<Description:> Get the output OBO filename

B<Parameters:> none

B<Returns:> $filename - scalar/string name of the output file

=cut

sub getFilename {

    my ($self) = shift;

    if (exists $self->{'_filename'}){
	return $self->{'_filename'};
    }
    else {
	$logger->warn("output filename does not exist");
    }
}

=item $obj->setDefaultNamespace()

B<Description:> Set the default-namespace for this ontology

B<Parameters:> default-namespace - scalar

B<Returns:>  none

=cut

sub setDefaultNamespace {

    my ($self) = shift;
    my ($defaultNamespace) = @_;

    if (!defined($defaultNamespace)){
	$logger->logdie("default-namespace was not defined");
    }

    $self->{'_headers'}->{'default-namespace'} = $defaultNamespace;
}

=item $obj->getDefaultNamespace()

B<Description:> Retrieve the default-namespace for this ontology

B<Parameters:> none

B<Returns:> default-namespace - scalar

=cut

sub getDefaultNamespace {

    my ($self) = shift;


    if (exists $self->{'_headers'}->{'default-namespace'}){
	return $self->{'_headers'}->{'default-namespace'};
    }
    elsif (exists $self->{'_headers'}->{'default_namespace'}){
	return $self->{'_headers'}->{'default_namespace'};
    }
    else {
	$logger->warn("default-namespace is not defined");
	return undef;
    }
}

=item $obj->doesHaveDefaultNamespace()

B<Description:> Determine whether a value for the default-namespace header attribute exists

B<Parameters:> none

B<Returns:> 

0 - false (scalar)
1 - true  (scalar)

=cut

sub doesHaveDefaultNamespace {

    my ($self) = shift;


    if (exists $self->{'_headers'}->{'default-namespace'}){
	return 1;
    }
    elsif (exists $self->{'_headers'}->{'default_namespace'}){
	return 1;
    }
    else {
	$logger->warn("default_namespace is not defined");
	return 0;
    }
}

=item $obj->setFormatVersion()

B<Description:> Set the format-version for this ontology

B<Parameters:> format-version - scalar

B<Returns:>  none

=cut

sub setFormatVersion {

    my ($self) = shift;
    my ($formatVersion) = @_;

    if (!defined($formatVersion)){
	$logger->logdie("format-version was not defined");
    }

    $self->{'_headers'}->{'format-version'} = $formatVersion;
}

=item $obj->getFormatVersion()

B<Description:> Retrieve the format-version for this ontology

B<Parameters:> none

B<Returns:> format-version - scalar

=cut

sub getFormatVersion {

    my ($self) = shift;


    if (exists $self->{'_headers'}->{'format-version'}){
	return $self->{'_headers'}->{'format-version'};
    }
    elsif (exists $self->{'_headers'}->{'format_version'}){
	return $self->{'_headers'}->{'format_version'};
    }
    else {
	$logger->warn("format-version is not defined");
	return undef;
    }
}

=item $obj->doesHaveFormatVersion()

B<Description:> Determine whether a value for the format-version header attribute exists

B<Parameters:> none

B<Returns:> 

0 - false (scalar)
1 - true  (scalar)

=cut

sub doesHaveFormatVersion {

    my ($self) = shift;


    if (exists $self->{'_headers'}->{'format-version'}){
	return 1;
    }
    elsif (exists $self->{'_headers'}->{'format_version'}){
	return 1;
    }
    else {
	$logger->warn("format-version is not defined");
	return 0;
    }
}

=item $obj->setDataVersion()

B<Description:> Set the data-version for this ontology

B<Parameters:> data-version - scalar

B<Returns:>  none

=cut

sub setDataVersion {

    my ($self) = shift;
    my ($dataVersion) = @_;

    if (!defined($dataVersion)){
	$logger->logdie("data-version was not defined");
    }

    $self->{'_headers'}->{'data-version'} = $dataVersion;
}

=item $obj->getDataVersion()

B<Description:> Retrieve the data-version for this ontology

B<Parameters:> none

B<Returns:> data-version - scalar

=cut

sub getDataVersion {

    my ($self) = shift;


    if (exists $self->{'_headers'}->{'data-version'}){
	return $self->{'_headers'}->{'data-version'};
    }
    elsif (exists $self->{'_headers'}->{'data_version'}){
	return $self->{'_headers'}->{'data_version'};
    }
    else {
	$logger->warn("data-version is not defined");
	return undef;
    }
}

=item $obj->doesHaveDataVersion()

B<Description:> Determine whether a value for the data-version header attribute exists

B<Parameters:> none

B<Returns:> 

0 - false (scalar)
1 - true  (scalar)

=cut

sub doesHaveDataVersion {

    my ($self) = shift;


    if (exists $self->{'_headers'}->{'data-version'}){
	return 1;
    }
    elsif (exists $self->{'_headers'}->{'data_version'}){
	return 1;
    }
    else {
	$logger->warn("data-version is not defined");
	return 0;
    }
}



=item $obj->setVersion()

B<Description:> Set the version for this ontology

B<Parameters:> version - scalar

B<Returns:>  none

=cut

sub setVersion {

    my ($self) = shift;
    my ($version) = @_;

    if (!defined($version)){
	$logger->logdie("version was not defined");
    }

    $self->{'_headers'}->{'version'} = $version;
}

=item $obj->getVersion()

B<Description:> Retrieve the version for this ontology

B<Parameters:> none

B<Returns:> version - scalar

=cut

sub getVersion {

    my ($self) = shift;


    if (exists $self->{'_headers'}->{'version'}){
	return $self->{'_headers'}->{'version'};
    }
    else {
	$logger->warn("version is not defined");
	return undef;
    }
}

=item $obj->doesHaveVersion()

B<Description:> Determine whether a value for the version header attribute exists

B<Parameters:> none

B<Returns:> 

0 - false (scalar)
1 - true  (scalar)

=cut

sub doesHaveVersion {

    my ($self) = shift;
    
    
    if (exists $self->{'_headers'}->{'version'}){
	return 1;
    }
    else {
	$logger->warn("version is not defined");
	return 0;
    }
}

=item $obj->setDate()

B<Description:> Set the date for this ontology

B<Parameters:> date - scalar

B<Returns:>  none

=cut

sub setDate {

    my ($self) = shift;
    my ($date) = @_;
    
    if (!defined($date)){
	$logger->logdie("date was not defined");
    }
    
    $self->{'_headers'}->{'date'} = $date;
}

=item $obj->getDate()

B<Description:> Retrieve the date for this ontology

B<Parameters:> none

B<Returns:> date - scalar

=cut

sub getDate {

    my ($self) = shift;


    if (exists $self->{'_headers'}->{'date'}){
	return $self->{'_headers'}->{'date'};
    }
    else {
	$logger->warn("date is not defined");
	return undef;
    }
}

=item $obj->doesHaveDate()

B<Description:> Determine whether a value for the date header attribute exists

B<Parameters:> none

B<Returns:> 

0 - false (scalar)
1 - true  (scalar)

=cut

sub doesHaveDateBy {

    my ($self) = shift;


    if (exists $self->{'_headers'}->{'date'}){
	return 1;
    }
    else {
	$logger->warn("date is not defined");
	return 0;
    }
}

=item $obj->setSavedBy()

B<Description:> Set the saved-by for this ontology

B<Parameters:> savedBy - scalar

B<Returns:>  none

=cut

sub setSavedBy {

    my ($self) = shift;
    my ($savedBy) = @_;

    if (!defined($savedBy)){
	$logger->logdie("savedBy was not defined");
    }

    $self->{'_headers'}->{'saved-by'} = $savedBy;
}

=item $obj->getSavedBy()

B<Description:> Retrieve the saved-by for this ontology

B<Parameters:> none

B<Returns:> saved-by - scalar

=cut

sub getSavedBy {

    my ($self) = shift;


    if (exists $self->{'_headers'}->{'saved-by'}){
	return $self->{'_headers'}->{'saved-by'};
    }
    elsif (exists $self->{'_headers'}->{'saved_by'}){
	return $self->{'_headers'}->{'saved_by'};
    }
    else {
	$logger->warn("saved-by is not defined");
	return undef;
    }
}

=item $obj->doesHaveSavedBy()

B<Description:> Determine whether a value for the saved-by header attribute exists

B<Parameters:> none

B<Returns:> 

0 - false (scalar)
1 - true  (scalar)

=cut

sub doesHaveSavedBy {

    my ($self) = shift;


    if (exists $self->{'_headers'}->{'saved-by'}){
	return 1;
    }
    elsif (exists $self->{'_headers'}->{'saved_by'}){
	return 1;
    }
    else {
	$logger->warn("saved-by is not defined");
	return 0;
    }
}

=item $obj->setAutoGeneratedBy()

B<Description:> Set the auto-generated-by for this ontology

B<Parameters:> autoGeneratedBy - scalar

B<Returns:>  none

=cut

sub setAutoGeneratedBy {

    my ($self) = shift;
    my ($autoGeneratedBy) = @_;

    if (!defined($autoGeneratedBy)){
	$logger->logdie("auto-generated-by was not defined");
    }

    $self->{'_headers'}->{'auto-generated-by'} = $autoGeneratedBy;
}

=item $obj->getAutoGeneratedBy()

B<Description:> Retrieve the auto-generated-by for this ontology

B<Parameters:> none

B<Returns:> auto-generated-by - scalar

=cut

sub getAutoGeneratedBy {

    my ($self) = shift;


    if (exists $self->{'_headers'}->{'auto-generated-by'}){
	return $self->{'_headers'}->{'auto-generated-by'};
    }
    elsif (exists $self->{'_headers'}->{'auto_generated_by'}){
	return $self->{'_headers'}->{'auto_generated_by'};
    }
    else {
	$logger->warn("auto-generated-by is not defined");
	return undef;
    }
}

=item $obj->doesHaveAutoGeneratedBy()

B<Description:> Determine whether a value for the auto-generated-by header attribute exists

B<Parameters:> none

B<Returns:> 

0 - false (scalar)
1 - true  (scalar)

=cut

sub doesHaveAutoGeneratedBy {

    my ($self) = shift;


    if (exists $self->{'_headers'}->{'auto-generated-by'}){
	return 1;
    }
    elsif (exists $self->{'_headers'}->{'auto_generated_by'}){
	return 1;
    }
    else {
	$logger->warn("subsetdef is not defined");
	return 0;
    }
}

=item $obj->setSubsetdef()

B<Description:> Set the subsetdef for this ontology

B<Parameters:> subsetdef - scalar

B<Returns:>  none

=cut

sub setSubsetdef {

    my ($self) = shift;
    my ($subsetdef) = @_;

    if (!defined($subsetdef)){
	$logger->logdie("subsetdef was not defined");
    }

    $self->{'_headers'}->{'subsetdef'} = $subsetdef;
}

=item $obj->getSubsetdef()

B<Description:> Retrieve the subsetdef for this ontology

B<Parameters:> none

B<Returns:> subsetdef - scalar

=cut

sub getSubsetdef {

    my ($self) = shift;


    if (exists $self->{'_headers'}->{'subsetdef'}){
	return $self->{'_headers'}->{'subsetdef'};
    }
    else {
	$logger->warn("subsetdef is not defined");
	return undef;
    }
}

=item $obj->doesHaveSubsetdef()

B<Description:> Determine whether a value for the subsetdef header attribute exists

B<Parameters:> none

B<Returns:> 

0 - false (scalar)
1 - true  (scalar)

=cut

sub doesHaveSubsetdef {

    my ($self) = shift;


    if (exists $self->{'_headers'}->{'subsetdef'}){
	return 1;
    }
    else {
	$logger->warn("subsetdef is not defined");
	return 0;
    }
}


=item $obj->setImport()

B<Description:> Set the import for this ontology

B<Parameters:> import - scalar

B<Returns:>  none

=cut

sub setImport {

    my ($self) = shift;
    my ($import) = @_;

    if (!defined($import)){
	$logger->logdie("import was not defined");
    }

    $self->{'_headers'}->{'import'} = $import;
}

=item $obj->getImport()

B<Description:> Retrieve the import for this ontology

B<Parameters:> none

B<Returns:> import - scalar

=cut

sub getImport {

    my ($self) = shift;


    if (exists $self->{'_headers'}->{'import'}){
	return $self->{'_headers'}->{'import'};
    }
    else {
	$logger->warn("import is not defined");
	return undef;
    }
}

=item $obj->doesHaveImport()

B<Description:> Determine whether a value for the import header attribute exists

B<Parameters:> none

B<Returns:> 

0 - false (scalar)
1 - true  (scalar)

=cut

sub doesHaveImport {

    my ($self) = shift;


    if (exists $self->{'_headers'}->{'import'}){
	return 1;
    }
    else {
	$logger->warn("import is not defined");
	return 0;
    }
}


=item $obj->setSynonymtypedef()

B<Description:> Set the synonymtypedef for this ontology

B<Parameters:> synonymtypedef - scalar

B<Returns:>  none

=cut

sub setSynonymtypedef {

    my ($self) = shift;
    my ($synonymtypedef) = @_;

    if (!defined($synonymtypedef)){
	$logger->logdie("synonymtypedef was not defined");
    }

    $self->{'_headers'}->{'synonymtypedef'} = $synonymtypedef;
}

=item $obj->getSynonymtypedef()

B<Description:> Retrieve the synonymtypedef for this ontology

B<Parameters:> none

B<Returns:> synonymtypedef - scalar

=cut

sub getSynonymtypedef {

    my ($self) = shift;


    if (exists $self->{'_headers'}->{'synonymtypedef'}){
	return $self->{'_headers'}->{'synonymtypedef'};
    }
    else {
	$logger->warn("synonymtypedef is not defined");
	return undef;
    }
}

=item $obj->doesHaveSynonymtypedef()

B<Description:> Determine whether a value for the synonymtypedef header attribute exists

B<Parameters:> none

B<Returns:> 

0 - false (scalar)
1 - true  (scalar)

=cut

sub doesHaveSynonymtypedef {

    my ($self) = shift;


    if (exists $self->{'_headers'}->{'synonymtypedef'}){
	return 1;
    }
    else {
	$logger->warn("synonymtypedef is not defined");
	return 0;
    }
}


=item $obj->setIdspace()

B<Description:> Set the idspace for this ontology

B<Parameters:> idspace - scalar

B<Returns:>  none

=cut

sub setIdspace {

    my ($self) = shift;
    my ($idspace) = @_;

    if (!defined($idspace)){
	$logger->logdie("idspace was not defined");
    }

    $self->{'_headers'}->{'idspace'} = $idspace;
}

=item $obj->getIdspace()

B<Description:> Retrieve the idspace for this ontology

B<Parameters:> none

B<Returns:> idspace - scalar

=cut

sub getIdspace {

    my ($self) = shift;


    if (exists $self->{'_headers'}->{'idspace'}){
	return $self->{'_headers'}->{'idspace'};
    }
    else {
	$logger->warn("idspace is not defined");
	return undef;
    }
}

=item $obj->doesHaveIdspace()

B<Description:> Determine whether a value for the idspace header attribute exists

B<Parameters:> none

B<Returns:> 

0 - false (scalar)
1 - true  (scalar)

=cut

sub doesHaveIdspace {

    my ($self) = shift;


    if (exists $self->{'_headers'}->{'idspace'}){
	return 1;
    }
    else {
	$logger->warn("idspace is not defined");
	return 0;
    }
}

=item $obj->setDefaultRelationshipIdPrefix()

B<Description:> Set the default-relationship-id-prefix for this ontology

B<Parameters:> defaultRelationshipIdPrefix - scalar

B<Returns:>  none

=cut

sub setDefaultRelationshipIdPrefix {

    my ($self) = shift;
    my ($defaultRelationshipIdPrefix) = @_;

    if (!defined($defaultRelationshipIdPrefix)){
	$logger->logdie("default-relationship-id-prefix was not defined");
    }

    $self->{'_headers'}->{'default-relationship-id-prefix'} = $defaultRelationshipIdPrefix;
}

=item $obj->getDefaultRelationshipIdPrefix()

B<Description:> Retrieve the default-relationship-id-prefix for this ontology

B<Parameters:> none

B<Returns:> default-relationship-id-prefix - scalar

=cut

sub getDefaultRelationshipIdPrefix {

    my ($self) = shift;


    if (exists $self->{'_headers'}->{'default-relationship-id-prefix'}){
	return $self->{'_headers'}->{'default-relationship-id-prefix'};
    }
    elsif (exists $self->{'_headers'}->{'default_relationship_id_prefix'}){
	return $self->{'_headers'}->{'default_relationship_id_prefix'};
    }
    else {
	$logger->warn("default-relationship-id-prefix is not defined");
	return undef;
    }
}

=item $obj->doesHaveDefaultRelationshipIdPrefix()

B<Description:> Determine whether a value for the default-relationship-id-prefix header attribute exists

B<Parameters:> none

B<Returns:> 

0 - false (scalar)
1 - true  (scalar)

=cut

sub doesHaveDefaultRelationshipIdPrefix {

    my ($self) = shift;


    if (exists $self->{'_headers'}->{'default-relationship-id-prefix'}){
	return 1;
    }
    elsif (exists $self->{'_headers'}->{'default-relationship-id-prefix'}){
	return 1;
    }
    else {
	$logger->warn("default-relationship-id-prefix is not defined");
	return 0;
    }
}

=item $obj->setIdMapping()

B<Description:> Set the id-mapping for this ontology

B<Parameters:> idMapping - scalar

B<Returns:>  none

=cut

sub setIdMapping {

    my ($self) = shift;
    my ($idMapping) = @_;

    if (!defined($idMapping)){
	$logger->logdie("id-mapping was not defined");
    }

    $self->{'_headers'}->{'id-mapping'} = $idMapping;
}

=item $obj->getIdMapping()

B<Description:> Retrieve the id-mapping for this ontology

B<Parameters:> none

B<Returns:> id-mapping - scalar

=cut

sub getIdMapping {

    my ($self) = shift;


    if (exists $self->{'_headers'}->{'id-mapping'}){
	return $self->{'_headers'}->{'id-mapping'};
    }
    elsif (exists $self->{'_headers'}->{'id_mapping'}){
	return $self->{'_headers'}->{'id_mapping'};
    }
    else {
	$logger->warn("id-mapping is not defined");
	return undef;
    }
}

=item $obj->doesHaveIdMapping()

B<Description:> Determine whether a value for the id-mapping header attribute exists

B<Parameters:> none

B<Returns:> 

0 - false (scalar)
1 - true  (scalar)

=cut

sub doesHaveIdMapping {

    my ($self) = shift;


    if (exists $self->{'_headers'}->{'id-mapping'}){
	return 1;
    }
    elsif (exists $self->{'_headers'}->{'id_mapping'}){
	return 1;
    }
    else {
	$logger->warn("id-mapping is not defined");
	return 0;
    }
}


=item $obj->setRemark()

B<Description:> Set the remark for this ontology

B<Parameters:> remark - scalar

B<Returns:>  none

=cut

sub setRemark {

    my ($self) = shift;
    my ($remark) = @_;

    if (!defined($remark)){
	$logger->logdie("remark was not defined");
    }

    $self->{'_headers'}->{'remark'} = $remark;
}

=item $obj->getRemark()

B<Description:> Retrieve the remark for this ontology

B<Parameters:> none

B<Returns:> remark - scalar

=cut

sub getRemark {

    my ($self) = shift;


    if (exists $self->{'_headers'}->{'remark'}){
	return $self->{'_headers'}->{'remark'};
    }
    else {
	$logger->warn("remark is not defined");
	return undef;
    }
}

=item $obj->doesHaveRemark()

B<Description:> Determine whether a value for the remark header attribute exists

B<Parameters:> none

B<Returns:> 

0 - false (scalar)
1 - true  (scalar)

=cut

sub doesHaveRemark {

    my ($self) = shift;


    if (exists $self->{'_headers'}->{'remark'}){
	return 1;
    }
    else {
	return 0;
    }
}


=item $obj->getHeaderByName($name)

B<Description:> Retrieves the header value by header name

B<Parameters:> $name - scalar/string name of the header

B<Returns:> $value - scalar/string

=cut

sub getHeaderByName {

    my ($self) = shift;
    my ($name) = @_;

    if (!defined($name)){
	$logger->logdie("name was not defined");
    }

    if (exists $self->{'_headers'}->{$name}){
	return $self->{'_headers'}->{$name};
    }
    else {
	$logger->warn("header with name '$name' does not exist");
    }
}

=item $obj->doesValueExistForHeader($name)

B<Description:> Determines whether a value exists for the named header attribute

B<Parameters:> $name - scalar/string name of the header

B<Returns:> 

 0 - false (scalar)
 1 - true  (scalar)

=cut

sub doesValueExistForHeader {

    my ($self) = shift;
    my ($name) = @_;

    if (!defined($name)){
	$logger->logdie("name was not defined");
    }

    if (exists $self->{'_headers'}->{$name}){
	return 1;
    }
    return 0;
}

=item $obj->createAndAddRecord( id=> $id, name => $name, def => $def, is_obsolete => $is_obsolete, type => $type)

B<Description:> Creates a new OBOTerm or OBOTypedefobject and stores reference to that object

B<Parameters:> 

$id          - scalar id
$name        - scalar name
$def         - scalar def
$is_obsolete - scalar is_obsolete
$type        - scalar type (term or typedef)

B<Returns:> reference to the OBOTerm or OBOTypedef

=cut

sub createAndAddRecord {

    my ($self) = shift;
    my (%args) = @_;

    if ((! exists $args{'type'}) || ($args{'type'} eq 'term')){
	return $self->createAndAddTerm(@_);
    }
    else {
	return $self->createAndAddTypedef(@_);
    }
}

sub createAndAddTerm {

    my ($self) = shift;
    my (%args) = @_;

    my $oboTerm = new OBO::OBOTerm($args{'id'}, $args{'name'}, $args{'def'}, $args{'is_obsolete'});
    if (!defined($oboTerm)){
	$logger->logdie("Could not instantiate OBO::OBOTerm for id '$args{'id'}' ".
			"name '$args{'name'}' ".
			"def '$args{'def'}' ".
			"is_obsolete '$args{'is_obsolete'}'");
    }
    
    ## Store reference to this new OBOTerm
    $self->{'_terms'}->{$args{'id'}} = $oboTerm;
    
    ## Increment the class member
    $oboTermCtr++;

    return $oboTerm;
}


sub createAndAddTypedef {

    my ($self) = shift;
    my (%args) = @_;

    my $oboTypedef = new OBO::OBOTypedef($args{'id'}, $args{'name'}, $args{'def'}, $args{'is_obsolete'});
    if (!defined($oboTypedef)){
	$logger->logdie("Could not instantiate OBO::OBOTypedef for id '$args{'id'}' name '$args{'name'}' def '$args{'def'}' is_obsolete '$args{'is_obsolete'}'");
    }	
    ## Store reference to this new OBOTypedef
    $self->{'_typedefs'}->{$args{'id'}} = $oboTypedef;

    ## Increment the class member
    $oboTypedef++;

    return $oboTypedef;
}



=item $obj->hasRecordWithId($id)

B<Description:> Check whether builder has reference to some OBOTerm or OBOTypedef with $id

B<Parameters:> 

$id - scalar/string corresponding to the id tag of the OBOTerm or OBOTypedef

B<Returns:> 

0 - scalar false
1 - scalar true

=cut

sub hasRecordWithId {
    my ($self) = shift;
    my ($id) = @_;

    if (!defined($id)){
	$logger->logdie("id was not defined");
    }

    if (exists $self->{'_terms'}->{$id}){
	return 1;
    }
    elsif (exists $self->{'_typedefs'}->{$id}){
	return 1;
    }
    else {
	return 0;
    }
}


=item $obj->getRecordById($id)

B<Description:> Return reference to OBOTerm or OBOTypedef with id $id

B<Parameters:> 

$id - scalar/string value that corresponds to the ID attribute of the OBOTerm or OBOTypedef

B<Returns:> reference to the OBOTerm or OBOTypedefwith id $id

=cut

sub getRecordById {

    my ($self) = shift;
    my ($id) = @_;

    if (!defined($id)){
	$logger->logdie("id was not defined");
    }
    
    if ( exists $self->{'_terms'}->{$id} ){
	return $self->{'_terms'}->{$id};
    }
    elsif ( exists $self->{'_typedefs'}->{$id} ){
	return $self->{'_typedefs'}->{$id};
    }
    else {
	$logger->warn("Neither OBOTerm nor OBOTypedef exists with id '$id'");
	return undef;
    }
}



=item $obj->getRecordNameById($id)

B<Description:> Return name attribute of OBOTerm or OBOTypedef with id $id

B<Parameters:> 

$id - scalar/string value that corresponds to the ID attribute of the OBOTerm or OBOTypedef

B<Returns:> name attribute of OBOTerm or OBOTypedefwith id $id

=cut

sub getRecordNameById {

    my ($self) = shift;
    my ($id) = @_;

    if (!defined($id)){
	$logger->logdie("id was not defined");
    }
    
    if ( exists $self->{'_terms'}->{$id} ){
	return $self->{'_terms'}->{$id}->getName();
    }
    elsif ( exists $self->{'_typedefs'}->{$id} ){
	return $self->{'_typedefs'}->{$id}->getName();
    }
    else {
	$logger->warn("Neither OBOTerm nor OBOTypedef exists with id '$id'");
	return undef;
    }
}

=item $obj->addTermById($id, $oboTerm)

B<Description:> Add reference to OBOTerm with id $id

B<Parameters:> 

$id       - scalar/string value that corresponds to the ID attribute of the OBOTerm
$oboTerm  - OBOTerm

B<Returns:> none

=cut

sub addTermById {

    my ($self) = shift;
    my ($id, $oboTerm) = @_;

    if (!defined($id)){
	$logger->logdie("id was not defined");
    }
    if (!defined($oboTerm)){
	$logger->logdie("oboTerm was not defined");
    }
    
    ## Ensure that no name-namespace tuples are repeated
    $self->checkForDuplicates($oboTerm);

    $self->{'_terms'}->{$id} = $oboTerm;

    ## The OBOTerm will maintain a reference to the OBOBuilder
    $oboTerm->{'_builder'} = $self;

    $oboTermCtr++;
}

=item $obj->addTypedefById($id, $oboTypedef)

B<Description:> Add reference to Typedef OBOTypedef with id $id

B<Parameters:> 

$id         - scalar/string value that corresponds to the ID attribute of the OBOTypedef
$oboTypedef - OBOTypedef

B<Returns:> none

=cut

sub addTypedefById {

    my ($self) = shift;
    my ($id, $oboTypedef) = @_;

    if (!defined($id)){
	$logger->logdie("id was not defined");
    }
    if (!defined($oboTypedef)){
	$logger->logdie("oboTypedef was not defined");
    }
    
    ## Ensure that no name-namespace tuples are repeated
    $self->checkForDuplicates($oboTypedef);

    $self->{'_typedefs'}->{$id} = $oboTypedef;

    $oboTypedefCtr++;
}

=item $obj->checkForDuplicates($oboTerm)

B<Description:> Check whether the name-namespace tuple has been repeated

B<Parameters:> OBOTerm or OBOTypedef

B<Returns:> none

=cut

sub checkForDuplicates {

    my $self = shift;
    my ($oboRecord) = @_;

    my $name = $oboRecord->getName();

    my $namespace;

    if ($oboRecord->hasNamespace()){
	$namespace = $oboRecord->getNamespace();
    }
    else {
	$namespace = $self->getDefaultNamespace();
    }

    if (exists $self->{'_name_namespace_tuples'}->{$name,$namespace}){
	$logger->logdie("name '$name' namespace '$namespace' already exist");
    }
}


=item $obj->removeRecordById($id)

B<Description:> Remove all references to the OBOTerm or OBOTypedef with id $id

B<Parameters:> $id (scalar)

B<Returns:> None

=cut

sub removeRecordById {

    my ($self) = shift;
    my ($id) = @_;

    if (exists $self->{'_terms'}->{$id}){
	delete $self->{'_terms'}->{$id};

	## Decrement the class member
	$oboTermCtr--;

	## Keep track of which id values for OBOTerms that have been removed.
	$self->{'_removed_terms'}->{$id}++;
    }
    elsif (exists $self->{'_typedefs'}->{$id}){
	delete $self->{'_typedefs'}->{$id};

	## Decrement the class member
	$oboTypedefCtr--;

	## Keep track of which id values for OBOTypedefs that have been removed.
	$self->{'_removed_typedefs'}->{$id}++;
    }
    else {
	$logger->warn("Neither OBOTerm nor OBOTypedef exists with id '$id'");
    }
}

sub writeFile {

    my ($self) = shift;

    if (exists $self->{'_filename'} ){
	my $fh;
	open ($fh, ">$self->{'_filename'}") || $logger->logdie("Could not open file '$self->{'_filename'}' in write mode: $!");
	$self->writeHeaders($fh);
	print $fh "\n\n";
	$self->writeTerms($fh);
	print $fh "\n";
	$self->writeTypedefs($fh);
	close $fh;	
    }
    else {
	$logger->logdie("output filename not defined");
    }
}

=item $obj->writeHeaders()

B<Description:> Writes the headers to the OBO file

B<Parameters:> 

$fh - file handle

B<Returns:> None

=cut

sub writeHeaders {

    my ($self) = shift;
    my ($fh) = @_;

    if (!defined($fh)){
	$logger->logdie("file handle was not defined");
    }

    if (exists $self->{'_headers'}){

	foreach my $header (@HEADER_ORDER){

	    if (exists $self->{'_headers'}->{$header}){

		print $fh "$header: $self->{'_headers'}->{$header}\n";
	    }
	    else {
		if (exists $REQUIRED_HEADERS->{$header}){
		    $logger->logdie("Required header '$header' was not defined");
		}
	    }
	}
    }
    else {
	$logger->logdie("_headers was not defined!");
    }
}

=item $obj->writeTerms()

B<Description:> Writes the OBOTerms to the OBO file

B<Parameters:> 

$fh - file handle

B<Returns:> None

=cut

sub writeTerms {

    my ($self) = shift;
    my ($fh) = @_;

    if (!defined($fh)){
	$logger->logdie("file handle was not defined");
    }
    
    foreach my $id ( sort keys %{$self->{'_terms'}} ) {
	$self->{'_terms'}->{$id}->writeRecord($fh);
	print $fh "\n";
    }
}


=item $obj->writeTypedefs()

B<Description:> Writes the OBOTypedefs to the OBO file

B<Parameters:> 

$fh - file handle

B<Returns:> None

=cut

sub writeTypedefs {

    my ($self) = shift;
    my ($fh) = @_;

    if (!defined($fh)){
	$logger->logdie("file handle was not defined");
    }
    
    foreach my $id (keys %{$self->{'_typedefs'}} ){
	$self->{'_typedefs'}->{$id}->writeRecord($fh);
	print $fh "\n";
    }
}

=item $obj->nextTerm()

B<Description:> Iteratively returns reference to each OBOTerm

B<Parameters:> None

B<Returns:> reference to OBOTerm

=cut

sub nextTerm {

    my ($self) = shift;

    if ($oboTermCtr > 0 ){

	if ($sorted == 0 ) {
	    if ($logger->is_debug()){
		$logger->debug("The OBOTerms have not yet been sorted by id");
	    }
	    ## If the OBOTerms aren't already sorted, do so now.
	    foreach my $id (sort keys %{$self->{'_terms'}} ) {
		push ( @{$self->{'_sorted_terms'}}, $id );	
	    }
	    
	    ## The list has been sorted once, let's not do it again.
	    $sorted = 1;
	}
	
	if ( $recordIndex < $oboTermCtr ){
	    my $id = $self->{'_sorted_terms'}->[$recordIndex++];
	    if (defined($id)){
		return $self->getRecordById($id);
	    }
	    else {
		$logger->logdie("id was not defined for recordIndex '$recordIndex'");
	    }
	}
    }
    else {
	return $self->nextTypedef();
    }

    return undef;
}

=item $obj->nextTypedef()

B<Description:> Iteratively returns reference to each OBOTypedef

B<Parameters:> None

B<Returns:> reference to OBOTypedef

=cut

sub nextTypedef {

    my ($self) = shift;

    if ($oboTypedefCtr > 0 ){

	if ($typedefs_sorted == 0 ) {
	    if ($logger->is_debug()){
		$logger->debug("The OBOTypedefs have not yet been sorted by id");
	    }
	    ## If the OBOTypedefs aren't already sorted, do so now.
	    foreach my $id (sort keys %{$self->{'_typedefs'}} ) {
		push ( @{$self->{'_sorted_typedefs'}}, $id );	
	    }
	    
	    ## The list has been sorted once, let's not do it again.
	    $typedefs_sorted = 1;
	}
	
	if ( $typedefRecordIndex < $oboTypedefCtr ){
	    my $id = $self->{'_sorted_typedefs'}->[$typedefRecordIndex++];
	    if (defined($id)){
		return $self->getRecordById($id);
	    }
	    else {
		$logger->logdie("id was not defined for recordIndex '$typedefRecordIndex'");
	    }
	}
    }
    
    return undef;
}

=item $obj->resetTermIndex()

B<Description:> Reset the OBOTerm index

B<Parameters:> None

B<Returns:> None

=cut

sub resetTermIndex {

    my ($self) = shift;
    $recordIndex=0;
}


=item $obj->getTermCount()

B<Description:> Get the number of Term records

B<Parameters:> none

B<Returns:> $count - scalar/string

=cut

sub getTermCount {

    my ($self) = shift;

    return $oboTermCtr;
}


1;

