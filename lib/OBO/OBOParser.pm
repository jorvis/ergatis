package OBO::OBOParser;
=head1 NAME

OBO::OBOParser.pm 

A class to facilitate the parsing of OBO files.

=head1 VERSION

1.0

=head1 SYNOPSIS

 use OBOParser;
 use OBOBuilder;
 my $parser = new OBO::OBOParser( filename => '/tmp/myfile.obo');
 my $builder = $parser->parseFile();

=head1 AUTHOR

Jay Sundaram
sundaram@jcvi.org

=head1 METHODS


=over 4

=cut

use strict;
use OBO::OBOBuilder;
use OBO::OBORecord;
use OBO::Logger;

my $logger = OBO::Logger::get_logger("Logger::OBO");

## Legal IDs and Special Identifiers
## Any string is a legal id, as long as it is not one of
## the built in identifiers. 
## Four of these are defined by the OBO spec.
## For the key-value pairs, the value has no meaning.
my $specialIdentifiers = { 'OBO:TYPE' => 1,
			   'OBO:TERM' => 1,
			   'OBO:TERM_OR_TYPE' => 1,
			   'OBO:INSTANCE' => 1 
		       };

## The OBO Flat File Format Specification version 1.2 (OFFFSV1.2) at
## http://www.geneontology.org/GO.format.obo-1_2.shtml
## states the rules for the OBO tags.
##
## In the comments that follow, the value in the key-value pair,
## descriptions have the following meaning:
##
## "required, processed" means that the particular tag is required 
## by OFFFSV1.2 and is processed by obov1p2tochado.
##
## "optional, processed" means that the particular tag is regarded
## as being optional by OFFFSV1.2 and all instances will be
## processed by obov1p2tochado.
##
## "optional, ignored" means that the particular tag is regarded
## as being optional by OFFFSV1.2 and no instances will be 
## processed by obov1p2tochado.

## Header Tags
##
## Cardinality rules: 
## Not explicitly stated in the OBO Flat File Format Specification version 1.2 (OFFFSV1.2).
## obov1p2tochado.pl will only accept one instance of each tag.
##
## For the key-value pairs, the value means:
## 1 : required, processed
## 2 : optional, processed
## 3 : optional, ignored
my $headerTags = { 'format-version' => 1,
		   'data-version' => 3,
		   'version' => 3,
		   'date' => 3,
		   'saved-by' => 3,
		   'auto-generated-by' => 3,
		   'subsetdef' => 3,
		   'import' => 3,
		   'synonymtypedef' => 3,
		   'idspace' => 3,
		   'default-relationship-id-prefix' => 2,
		   'id-mapping' => 3,
		   'remark'  => 3,
		   'default-namespace' => 2
	       };    


## There are only three stanza types
##
## For the key-value pairs, the value means:
## 1 : required, processed
## 2 : optional, processed
## 3 : optional, ignored
my $stanzaTypes = { 'Typedef' => 2,
		    'Term' => 1,
		    'Instance' => 3
		};



## Term Stanza Tags
##
## For the key-value pairs, the value means:
## 1 : required, processed
## 2 : optional, processed
## 3 : optional, ignored
my $termStanzaTags = { 'id' => 1,
		       'name' => 1,
		       'is_anonymous' => 3,
		       'alt_id' => 2,
		       'def' => 2,
		       'comment' => 2,
		       'subset' => 3,
		       'synonym' => 2,
		       'exact_synonym' => 2,
		       'related_synonym' => 2,
		       'narrow_synonym' => 2,
		       'broad_synonym' => 2,
		       'xref' => 2,
		       'is_a' => 2,
		       'intersection_of' => 3,
		       'union_of' => 3,
		       'disjoint_from' => 3,
		       'relationship' => 2,
		       'is_obsolete' => 2,
		       'replaced_by' => 3,
		       'consider' => 3,
		       'builtin' => 3,
		       'namespace' => 2,
		       'use_term' => 3,    ## OBO1.0, deprecated, ignored
		       'xref_analog' => 2, ## OBO1.0, deprecated but support here
		       'xref_unk' => 2     ## OBO1.0, deprecated but support here
		   };

## Typedef Stanza Tags
##
## For the key-value pairs, the value means:
## 1 : required, processed 
## 2 : optional, processed
## 3 : optional, ignored
my $typedefStanzaTags = { 'id' => 1,
			  'name' => 1,
			  'is_anonymous' => 3,
			  'alt_id' => 2,
			  'def' => 2,
			  'comment' => 2,
			  'subset' => 3,
			  'synonym' => 2,
			  'exact_synonym' => 2,
			  'related_synonym' => 2,
			  'narrow_synonym' => 2,
			  'broad_synonym' => 2,
			  'xref' => 2,
			  'is_a' => 2,
			  'relationship' => 2,
			  'is_obsolete' => 2,
			  'replaced_by' => 3,
			  'consider' => 3,
			  'builtin' => 3,     ## Typedef only tags follow
			  'domain' => 3,       
			  'range' => 3,
			  'inverse_of' => 3,
			  'transitive_over' => 3,
			  'is_cyclic' => 3,
			  'is_reflexive' => 3,
			  'is_symmetric' => 3,
			  'is_anti_symmetric' => 3,
			  'is_transitive' => 3,
			  'is_metadata_tag' => 3,
			  'use_term' => 3,    ## OBO1.0, deprecated, ignored
			  'xref_analog' => 2, ## OBO1.0, deprecated but support here
			  'xref_unk' => 2     ## OBO1.0, deprecated but support here
		      };


## Cardinality rules for the Typedef Stanza Tags
##
## For the key-value pairs, the value means:
## 1: Only 1
## 2: Only 2 (etc.)
## m: one-to-many
my $typedefStanzaTagCardinality = { 'id' => 1,
				    'name' => 1,
				    'is_anonymous' => 1,
				    'alt_id' => 'm',
				    'def' => 1,
				    'comment' => 1,
				    'subset' => 1,
				    'synonym' => 1,
				    'exact_synonym' => 1,
				    'related_synonym' => 1,
				    'narrow_synonym' => 1,
				    'broad_synonym' => 1,
				    'xref' => 'm',
				    'is_a' => 1,
				    'relationship' => 1,
				    'is_obsolete' => 1,
				    'replaced_by' => 1,
				    'consider=' => 1,
				    'builtin' => 1,
				    'domain' => 1,
				    'range' => 1,
				    'inverse_of' => 1,
				    'transitive_over' => 1,
				    'is_cyclic' => 1,
				    'is_reflexive' => 1,
				    'is_symmetric' => 1,
				    'is_anti_symmetric' => 1,
				    'is_transitive' => 1,
				    'is_metadata_tag' => 1,
				    'xref_analog' => 'm',
				};


## Cardinality rules for the Term Stanza Tags
##
##
## For the key-value pairs, the value means:
## 1: Only 1
## 2: Only 2 (etc.)
## m: one-to-many
my $termStanzaTagCardinality = { 'id' => 1,
				 'name' => 1,
				 'is_anonymous' => 1,
				 'alt_id' => 'm',
				 'def' => 1,
				 'comment' => 1,
				 'subset' => 1,
				 'synonym' => 'm',
				 'exact_synonym' => 'm',
				 'related_synonym' => 'm',
				 'narrow_synonym' => 'm',
				 'broad_synonym' => 'm',
				 'xref' => 'm',
				 'is_a' => 'm',
				 'intersection_of' => 'm',
				 'union_of' => 'm',
				 'disjoint_from' => 'm',
				 'relationship' => 'm',
				 'is_obsolete' => 1,
				 'replaced_by' => 1,
				 'consider=' => 1,
				 'builtin' => 1,
				 'namespace' => 1,
				 'xref_analog' => 'm'
			     };


## Instance Stanza Tags
##
## For the key-value pairs, the value means:
## 1 : required, processed
## 2 : optional, processed
## 3 : optional, ignored
my $instanceOfStanzaTags = { 'id' => 1,
			     'name' => 1,
			     'instance_of' => 1,
			     'property_value' => 3,
			     'is_anonymous' => 3,
			     'namespace' => 3,
			     'alt_id' => 3,
			     'comment' => 3,
			     'xref' => 3,
			     'synonym' => 3,
			     'is_obsolete' => 3,
			     'replaced_by' => 3,
			     'consider' => 3,
			     'xref_analog' => 3,
			 };



## If Term Stanza contains tag is_obsolete == true
## then the obsoleted term should not have the 
## following tags.
## For the key-value pairs, the value has no meaning.
my $ifIsObsoleteDisallowTags = { 'is_a' => 1,
				 'relationship' => 1,
				 'inverse_of' => 1,
				 'disjoint_from' => 1,
				 'union_of' => 1,
				 'intersection_of' => 1
			     };


# Keep count of all tags encountered
my $oboTagCounter = {};

# Keep track of all unique id values encountered
my $uniqueIdLookup = {};

## Keep track of unexpected tags
my $unexpectedHeaderTags = {};
my $unexpectedTermStanzaTags = {};
my $unexpectedTypedefStanzaTags = {};
my $unexpectedInstanceStanzaTags = {};


## Keep track of ignored tags
my $ignoredHeaderTags = {};
my $ignoredTermStanzaTags = {};
my $ignoredTypedefStanzaTags = {};
my $ignoredInstanceStanzaTags = {};

## Keep track of all encountered tags
my $headerTagsTracker = {};
my $termStanzaTagsTracker = {};
my $typedefStanzaTagsTracker = {};
my $instanceStanzaTagsTracker = {};

=item new()

B<Description:> Instantiate OBOParser object

B<Parameters:> 

 %args

B<Returns:> Returns a reference to OBOParser

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
}

=item DESTROY

B<Description:> OBOParser class destructor

B<Parameters:> None

B<Returns:> None

=cut

sub DESTROY  {
    my $self = shift;

}

=item $obj->setFilename($filename)

B<Description:> Set the input OBO filename

B<Parameters:> $filename - scalar/string name of the input file

B<Returns:> None

=cut

sub setFilename {

    my ($self) = shift;
    my ($filename) = @_;

    if (!defined($filename)){
	$logger->logdie("input filename was not defined");
    }

    $self->{'_filename'} = $filename;
}

=item $obj->getFilename($filename)

B<Description:> Get the input OBO filename

B<Parameters:> none

B<Returns:> $filename - scalar/string name of the input file

=cut

sub getFilename {

    my ($self) = shift;

    if (exists $self->{'_filename'}){
	return $self->{'_filename'};
    }
    else {
	$logger->warn("input filename does not exist");
    }
}

=item $obj->parseHeadersOnly()

B<Description:> Parses only the header section of the input OBO file and creates/returns an OBOBuilder object

B<Parameters:> None

B<Returns:> reference to OBOBuilder

=cut

sub parseHeadersOnly {

    my $self = shift;

    if (exists $self->{'_filename'}){

	my $oboBuilder = new OBO::OBOBuilder();

	return $oboBuilder;
    }
    else {
	$logger->logdie("filename was not defined");
    }
}


=item $obj->parseFile()

B<Description:> Parses the input OBO file and creates/returns an OBOBuilder object

B<Parameters:> None

B<Returns:> reference to OBOBuilder

=cut

sub parseFile {

    my $self = shift;

    if (exists $self->{'_filename'}){

	## Flag for tracking when any one of the three stanza types have been encountered.
	my $someStanzaEncountered = 0;

	# keep track of line numbers
	my $lineCounter=0;

	my $fileContents = $self->_getFileContents($self->{'_filename'});

	print "Processing '$self->{'_filename'}' contents\n";

	#-----------------------------------------------------------------
	# show_progress related data
	#
	#----------------------------------------------------------------
# 	my $row_count = 0;
# 	my $bars = 30;
# 	my $total_rows = scalar(@{$fileContents});
# 	my $counter = int(.01 * $total_rows);
# 	$counter = 1 if($counter ==0);

# 	my $progressReporter = new ProgressReporter( bars => 30,
# 						     total_rows => scalar(@{$fileContents}) );

	my $oboBuilder = new OBO::OBOBuilder();

       	my $oboTerm; ## For reference to OBO::OBOTerm object

	my $oboTypedef; ## For reference to OBO::OBOTypedef object

	# The current unique id variable is defined outside the loop
	my $id;

	## Keep track of the current type of stanza that is being processed
	my $stanzaType;

	foreach my $line (@{$fileContents}){

	    $lineCounter++;
#	    $progressReporter->increment();

	    if ($line =~ /^\s*$/){
		next; # skip blank lines
	    }

	    if ($line =~ /^!/){
		next; # skip commented lines
	    }

# 	    $row_count++;

# 	    if ($logger->info("parsing")) {
# 		$prism->show_progress("Parsing OBO file $row_count/$total_rows",
# 				      $counter,$row_count,$bars,$total_rows) 
# 	    }

	    if ($line =~ /^\[(.*)\]/){
		## We've encountered the beginning of the next stanza

		$stanzaType = $1;
		
		$someStanzaEncountered = 1;

		if ( exists $stanzaTypes->{$stanzaType} ){
		    
		    ## Figure out what kind of stanza we're processing.
		    ## This will have bearing on how the tags in this
		    ## stanza will be processed.
		    ##
		    if ($stanzaType eq 'Typedef'){
			if (defined($oboTypedef)){
			    ## All data related to the previously encountered Typedef
			    ## has been processed.  Add the OBO::OBOTypedef to the
			    ## OBO::OBOBuilder.
			    $oboBuilder->addTypedefById($id, $oboTypedef);
			}
			if (defined($oboTerm)){
			    ## All data related to the previously encountered Term
			    ## has been processed.  Add the OBO::OBOTerm to the
			    ## OBO::OBOBuilder.
			    $oboBuilder->addTermById($id, $oboTerm);
			}

			## Create a new OBO::OBOTypedef for the next Typedef
			$oboTypedef = new OBO::OBOTypedef('id', 'name');
		    }
		    elsif ($stanzaType eq 'Term'){
			if (defined($oboTerm)){
			    ## All data related to the previously encountered Term
			    ## has been processed.  Add the OBO::OBOTerm to the
			    ## OBO::OBOBuilder.
			    $oboBuilder->addTermById($id, $oboTerm);
			}

			## Create a new OBO::OBOTerm for the next Term
			$oboTerm = new OBO::OBOTerm( 'id', 'name');

		    }
		    elsif ($stanzaType eq 'Instance'){
			if ($logger->is_debug()){
			    $logger->debug("Instance stanza encountered");
			}
		    }
		    else {
			$logger->logdie("Unrecognized stanza type '$stanzaType'");
		    }
		}
		else {
		    $logger->logdie("Unexpected stanza type '$stanzaType' encountered at line '$lineCounter' in OBO file '$self->{'_filename'}'");
		}

		next; ## Go ahead and read the next line
	    }

	    if ($line =~ /^(\S+):(.+)$/){
		## Parse all tag:value pairs

		my $tag = $1;

		my $value = $2;

		if (!defined($value)){
		    $logger->logdie("value was not defined for tag '$tag' while processing id '$id'");
		}

		$value =~ s/^\s+//;
		$value =~ s/\s+$//;
		
		## Maintain a count of all tag types encountered.
		$oboTagCounter->{$tag}++;

		if ($someStanzaEncountered == 0) {
		    ## Have not yet encountered a stanza therefore still processing header section
		    $self->_processHeaderData($oboBuilder, $tag, $value, $lineCounter);
		}
		elsif ( $stanzaType eq 'Term' ){
		    ## Process the Term stanza
		    $id = $self->_processTermStanza($oboTerm, $tag, $value, $id, $lineCounter);
		}
		elsif ( $stanzaType eq 'Typedef' ){
		    ## Process this Typedef stanza
		    $id = $self->_processTypedefStanza($oboTypedef, $tag, $value, $id, $lineCounter);
		}
 		elsif ( $stanzaType eq 'Instance' ) {
		    ## Process this Instance stanza

		    ## Keep track of all Instance stanza tags encountered
		    $instanceStanzaTagsTracker->{$tag}++;
		    
		    $logger->warn("$0 currently does not process Instance Stanza tags");
		}
		else {
		    $logger->logdie("tag '$tag' value '$value' ".
				    "someStanzaEncountered '$someStanzaEncountered' ".
				    "stanzaType '$stanzaType' ".
				    "at line '$lineCounter' in OBO file '$self->{'_filename'}'");
		}
	    }
	    else {
		$logger->logdie("Unexpected content at line '$lineCounter': $line");
	    }
	}

	if ($someStanzaEncountered == 1){
	    ## The last record in the file

	    if ($stanzaType eq 'Typedef'){
		if (defined($oboTypedef)){
		    ## All data related to the previously encountered Typedef
		    ## has been processed.  Add the OBO::OBOTypedef to the
		    ## OBO::OBOBuilder.
		    print STDERR $id . "\n";
		    $oboBuilder->addTypedefById($id, $oboTypedef);
		}
		
		## Create a new OBO::OBOTypedef for the next Typedef
		$oboTypedef = new OBO::OBOTypedef('id', 'name');
	    }
	    elsif ($stanzaType eq 'Term'){
		if (defined($oboTerm)){
		    ## All data related to the previously encountered Term
		    ## has been processed.  Add the OBO::OBOTerm to the
		    ## OBO::OBOBuilder.
		    $oboBuilder->addTermById($id, $oboTerm);
		}
		
		## Create a new OBO::OBOTerm for the next Term
		$oboTerm = new OBO::OBOTerm( 'id', 'name');
		
	    }
	    elsif ($stanzaType eq 'Instance'){
		if ($logger->is_debug()){
		    $logger->debug("Instance stanza encountered");
		}
	    }
	    else {
		$logger->logdie("Unrecognized stanza type '$stanzaType'");
	    }
	}

	$logger->info("Encountered '$oboTagCounter->{'name'}' term name values in OBO file '$self->{'_filename'}'");

	$self->{'_file_parsed'} = 1;

	return $oboBuilder;
    }
    else {
	$logger->logdie("filename was not defined");
    }


}

=item postParsingReport()

B<Description:> Report on the tags encounter during the file parse

B<Parameters:> none

B<Returns:> none

=cut

sub postParsingReport {

    my $self = shift;

    if ($self->{'_file_parsed'} == 1){

	## Reports for Header tags
	&_reportTagsEncountered($headerTagsTracker, $headerTags, "Header");
	&_reportUnexpectedTagsEncountered($unexpectedHeaderTags, "Header");
	&_reportIgnoredTagsEncountered($ignoredHeaderTags, "Header");

	## Reports for Term stanza tags
	&_reportTagsEncountered($termStanzaTagsTracker, $termStanzaTags, "Term");
	&_reportUnexpectedTagsEncountered($unexpectedTermStanzaTags, "Term");
	&_reportIgnoredTagsEncountered($ignoredTermStanzaTags, "Term");

	## Reports for Typedef stanza tags
	&_reportTagsEncountered($typedefStanzaTagsTracker, $typedefStanzaTags, "Typedef");
	&_reportUnexpectedTagsEncountered($unexpectedTypedefStanzaTags, "Typedef");
	&_reportIgnoredTagsEncountered($ignoredTypedefStanzaTags, "Typedef");

	## Reports for Instance tags
	&_reportTagsEncountered($instanceStanzaTagsTracker, $instanceOfStanzaTags,  "Instance");
    }
    else {
	$logger->warn("Can't report anything since the file '$self->{'_filename'}' has not yet been parsed!");
    }
}

=item _processHeaderData()

B<Description:> Process the header data

B<Parameters:> 

  $oboBuilder - OBOBuilder
  $tag        - scalar name of the header
  $value      - scalar value

B<Returns:> none

=cut

sub _processHeaderData {

    my $self = shift;
    my ($oboBuilder, $tag, $value, $lineCounter) = @_;

    ## Keep track of all header tags encountered
    $headerTagsTracker->{$tag}++;

    if ( exists $headerTags->{$tag} ) {
	## This is a recognized tag
	
	if ( $headerTags->{$tag} != 3 ){
	    ## This tag-value shall be processed

	    if ( $oboBuilder->doesHeaderExist($tag)){
		$logger->warn("Already added header with tag '$tag'");
	    }
	    else {
		$oboBuilder->addHeader($tag, $value);
	    }
	}
	else {
	    $ignoredHeaderTags->{$tag}++;
	}
    }
    else {
	## This tag is unrecognized
	$logger->warn("Unexpected Header tag '$tag' encountered at line '$lineCounter' in OBO file '$self->{'_filename'}'");
	$unexpectedHeaderTags->{$tag}++;
    }
}

=item _processTermStanza()

B<Description:> Process term stanza data

B<Parameters:> 

  $oboTerm - OBOTerm
  $tag     - scalar tag name
  $value   - scalar value
  $id      - scalar id
  $lineCounter - scalar line number in the OBO file

B<Returns:> $id - scalar

=cut

sub _processTermStanza {

    my $self = shift;
    my ($oboTerm, $tag, $value, $id, $lineCounter) = @_;

    ## Keep track of all Term stanza tags encountered
    $termStanzaTagsTracker->{$tag}++;

    if ( exists $termStanzaTags->{$tag} ){
	## This is a recognized tag
	if ( $termStanzaTags->{$tag} != 3){
	    ## This tag shall be processed
	    $id = $self->_storeTermTagValue($oboTerm, 
					    $tag,
					    $value, 
					    $id, 
					    $lineCounter);
	}
	else {
	    ## This tag was ignored
	    $ignoredTermStanzaTags->{$tag}++;
	}
    }
    else {
	## This tag is unrecognized
	$logger->warn("Unexpected Term stanza tag '$tag' encountered at line '$lineCounter'");
	$unexpectedTermStanzaTags->{$tag}++;
    }

    return $id;
}

=item processTypedefStanza()

B<Description:> Process typedef stanza data

B<Parameters:> 

  $oboTerm - OBOTypedef
  $tag     - scalar tag name
  $value   - scalar value
  $id      - scalar id
  $lineCounter - scalar line number in the OBO file

B<Returns:> $id - scalar

=cut

sub _processTypedefStanza {

    my $self = shift;
    my ($oboTypedef, $tag, $value, $id, $lineCounter) = @_;

    ## Keep track of all Typedef stanza tags encountered
    $typedefStanzaTagsTracker->{$tag}++;

    if ( exists $typedefStanzaTags->{$tag} ){
	## This is a recognized tag

	if ($typedefStanzaTags->{$tag} != 3 ){
	    ## This tag shall be processed
	    $id = $self->_storeTermTagValue($oboTypedef,
					    $tag,
					    $value,
					    $id,
					    $lineCounter);
	}
	else {
	    ## This tag was ignored
	    $ignoredTypedefStanzaTags->{$tag}++;
	}
    }
    else {
	## This tag is unrecognized
	$logger->warn("Unexpected Typedef stanza tag '$tag' encountered at line '$lineCounter'");
	$unexpectedTypedefStanzaTags->{$tag}++;
    }

    return $id;
}

=item $obj->_getFileContents()

B<Description:> Read in entire contents of the OBO file

B<Parameters:> $filename - scalar name of the OBO file

B<Returns:> Reference to array

=cut

sub _getFileContents {

    my $self = shift;
    my ($filename) = @_;

    if (!defined($filename)){
	$logger->logdie("filename was not defined");
    }

    if (!-e $filename){
	$logger->logdie("file '$filename' does not exist");
    }

    if (!-r $filename){
	$logger->logdie("file '$filename' does not have read permissions");
    }

    if (!-f $filename){
	$logger->logdie("file '$filename' is not a regular file");
    }

    if (!-s $filename){
	$logger->logdie("file '$filename' does not have any content");
    }


    print "Reading entire contents of file '$filename'\n";

    open (INFILE, "<$filename") or $logger->logdie("Could not open file '$filename' in read mode: $!");

    my @lines = <INFILE>;
    
    chomp @lines;

    return \@lines;
}

=item setTermStanzaTags()

B<Description:> TBA

B<Parameters:> 

  $termStanzaTags - reference to hash

B<Returns:> none

=cut

sub setTermStanzaTags {

    my ($termStanzaTags, $ignoreRelationships, $relationshipsOnly) = @_;

    if ($ignoreRelationships == 1){

	## These tags become optional, ignored
	$termStanzaTags->{'is_a'} = 3;
	$termStanzaTags->{'relationship'} = 3;

    }	
    if ($relationshipsOnly == 1){

	## All tags save id, is_a, relationship,
	## namespace become optional, ignored
	$termStanzaTags = { 'id' => 1,
			    'name' => 3,
			    'is_anonymous' => 3,
			    'alt_id' => 3,
			    'def' => 3,
			    'comment' => 3,
			    'subset' => 3,
			    'synonym' => 3,
			    'exact_synonym' => 3,
			    'related_synonym' => 3,
			    'narrow_synonym' => 3,
			    'broad_synonym' => 3,
			    'xref' => 3,
			    'is_a' => 2,
			    'intersection_of' => 3,
			    'union_of' => 3,
			    'disjoint_from' => 3,
			    'relationship' => 2,
			    'is_obsolete' => 3,
			    'replaced_by' => 3,
			    'consider=' => 3,
			    'builtin' => 3,
			    'namespace' => 2,
			    'xref' => 3
			   };
    }

}

=item _storeTermTagValue()

B<Description:> Extracts tag-value from line and adds to the OBOTerm or OBOTypedef

B<Parameters:> 

  $oboTerm - OBOTerm or OBOTypedef
  $tag     - scalar
  $value   - scalar
  $id      - scalar
  $lineCounter - scalar

B<Returns:> none

=cut

sub _storeTermTagValue {

    my $self = shift;
    my ($oboTerm, $tag, $value, $id, $lineCounter) = @_;

    if ($tag eq 'id'){

	if ( exists $specialIdentifiers->{$value} ){
	    $logger->logdie("Encountered reserved id '$value' at line '$lineCounter' in OBO file '$self->{'_filename'}'");
	}

	if ( exists $uniqueIdLookup->{$value} ){
	    $logger->logdie("Found duplicate id '$value' at line '$lineCounter' of OBO file '$self->{'_filename'}'");
	}

	$uniqueIdLookup->{$value}++;

	$oboTerm->setId($value);
	return $value;
    }
    else {
	if (!defined($id)){
	    $logger->logdie("id is not defined for tag '$tag' value '$value'");
	}

	if (!defined($value)){
	    $logger->logdie("value was not defined for tag '$tag' while processing id '$id'");
	}


	if ($tag eq 'name'){
	    $oboTerm->setName($value);
	}
	elsif ($tag eq 'namespace'){
	    $oboTerm->setNamespace($value);
	}
	elsif (($tag eq 'is_obsolete') && ($value eq 'true')){
	    $oboTerm->setIsObsolete(1);
	}
	elsif ($tag eq 'is_a'){
	    $oboTerm->addIsA(&_remove_trailing_comment($value));
	}
	elsif ($tag eq 'relationship'){
	    &_storeRelationship($oboTerm, $id, $value);
	}
	elsif ($tag eq 'def'){
	    $oboTerm->setDef(&_extract_value_from_double_quotes($value));
	}
	elsif ($tag eq 'comment'){
	    $oboTerm->setComment(&_extract_value_from_double_quotes($value));
	}
	elsif ($tag eq 'alt_id'){
	    $oboTerm->addAltId($value);
	}
	elsif ($tag eq 'xref_analog'){
	    $oboTerm->addXref($value);
	}
	elsif ($tag =~ /synonym/){
	    if ($value =~ /\"(.+)\"\s+(.+)/){
		my $synonymString = $1;
		my $optionalInfo = $2;

		$value = $synonymString;
		if ($optionalInfo =~ /EXACT/){
		    $tag = 'exact_synonym';
		}
		elsif ($optionalInfo =~ /NARROW/){
		    $tag = 'narrow_synonym';
		}
		elsif ($optionalInfo =~ /BROAD/){
		    $tag = 'broad_synonym';
		}
		elsif ($optionalInfo =~ /RELATED/){
		    $tag = 'related_synonym';
		}
	    }
	    else {
		$logger->logdie("Could not parse synonym type value '$value' with tag '$tag'");
	    }
	}
	else {
	    $logger->logdie("Encountered unexpected tag '$tag' with value '$value' while processing id '$id'");
	}
    }

    return $id;
}

=item _extract_value_from_double_quotes()

B<Description:> Extracts value from double quotes

B<Parameters:> $value - scalar

B<Returns:> $value - scalar

=cut

sub _extract_value_from_double_quotes {

    my ($value) = @_;
    
    if  (($value =~ /^\"/) && ($value =~ /\"$/)){

	if ($value =~ /\"(.+)\"/){
	    
	    my $newvalue = $1;
	    
	    # remove leading whitespaces
	    $newvalue =~ s/^\s*//;
	    
	    # remove trailing whitespaces
	    $newvalue =~ s/\s*$//;
	    
	    return $newvalue;
	}
	else {
	    $logger->logdie("Could not parse synonym '$value'");
	}
    }
    else {
#	print STDERR "my comment: $value\n";
	return $value;
    }
}

=item _remove_trailing_comment()

B<Description:> TBA

B<Parameters:> $value - scalar

B<Returns:> $value - scalar

=cut

sub _remove_trailing_comment {

    my ($value) = @_;

    if ($value =~ /^(\S+)\s*/){
	$value = $1;
    }
    else {
	$logger->logdie("Could not parse value '$value'");
    }
    
    return $value;
}

=item _storeRelationship()

B<Description:> Add relationship information to the OBOTerm or OBOTypedef

B<Parameters:> 

  $oboTerm - OBOTerm or OBOTypedef
  $id      - scalar
  $value   - scalar

B<Returns:> none

=cut

sub _storeRelationship {

    my ($oboTerm, $id, $value) = @_;

    # E.g.:
    # relationship: part_of SO:0000673 ! transcript

    if (!defined($id)){
	$logger->logdie("id was not defined");
    }

    if (!defined($value)){
	$logger->logdie("value was not defined");
    }

    my ($reltype, $val) = split(/\s+/, $value);

    if (!defined($reltype)){
	$logger->logdie("reltype not defined for value '$value' id '$id'");
    }

    if (!defined($val)){
	$logger->logdie("reltype not defined for value '$value' id '$id'");
    }
    
    $oboTerm->addRelationship($reltype, $val);
}

=item _reportTagsEncountered()

B<Description:> Report to log file a tally of all tags encountered

B<Parameters:> 

  $encounteredTagsLooup - reference to hash
  $tagsLookup           - reference to hash
  $type                 - scalar

B<Returns:> none

=cut

sub _reportTagsEncountered {

    my ($encounteredTagsLookup, $tagsLookup, $type) = @_;

    $logger->warn("Examining the '$type' tags");

    foreach my $tag (keys %{$tagsLookup}){
	if (exists $encounteredTagsLookup->{$tag}){
	    $logger->warn("Found '$encounteredTagsLookup->{$tag}' instances of '$tag'");
	}
	else {
	    $logger->warn("Found no occurences of '$tag'");
	}
    }

    foreach my $tag (keys %{$encounteredTagsLookup}){
	if (!exists $tagsLookup->{$tag}){
	    $logger->warn("Found '$encounteredTagsLookup->{$tag}' instances of '$tag'.  These tags are not part of the OBO 1.2 spec.");
	}
    }
}

=item _reportUnexpectedTagsEncountered()

B<Description:> Report to log file a tally of all unexpected tags encountered

B<Parameters:> 

  $unexpectedTagsLookup - reference to hash
  $type                 - scalar

B<Returns:> none

=cut

sub _reportUnexpectedTagsEncountered {

    my ($unexpectedTagsLookup, $type) = @_;

    $logger->warn("Examining unexpected '$type' tags that were encountered. Note that $0 will not process these tags.");

    foreach my $tag (keys %{$unexpectedTagsLookup}){
	$logger->warn("Unexpectedly encountered '$unexpectedTagsLookup->{$tag}' instances of tag '$tag'");
    }
}

=item _reportIgnoredTagsEncountered()

B<Description:> Report to log file a tally of all ignored tags encountered

B<Parameters:> 

  $unexpectedTagsLookup - reference to hash
  $type                 - scalar

B<Returns:> none

=cut

sub _reportIgnoredTagsEncountered {

    my ($ignoredTagsLookup, $type) = @_;

    $logger->warn("Examining ignored '$type' tags that were encountered. Note that $0 will not process these tags.");

    foreach my $tag (keys %{$ignoredTagsLookup}){
	$logger->warn("'$ignoredTagsLookup->{$tag}' instance(s) of tag '$tag' will be ignored");
    }
}


1;

