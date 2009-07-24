package GFF3::GFF3Builder;
=head1 NAME

GFF3::GFF3Builder.pm 

A class to facilitate the creation of GFF3 files.

=head1 VERSION

1.0

=head1 SYNOPSIS

 use GFF3Builder;
 my $obj = new GFF3::GFF3Builder( filename => '/tmp/myfile.gff3');
 $obj->addRecord($gff3Record);
 $obj->writeRecords();

=head1 AUTHOR

Jay Sundaram
sundaram@tigr.org

=head1 METHODS


=over 4

=cut

use strict;
use GFF3::GFF3Record;
use GFF3::Logger;

my $logger = GFF3::Logger::get_logger("Logger::GFF3");


## Class variable for keeping track of the number of GFF3 records
my $gff3RecordCtr = 0;

## Class variable for keeping track of the number of GFF3 records with FASTA sequence.
my $fastaRecordCtr = 0;

## Instance variable for supporting the nextRecord() method
my $recordIndex = 0;

## Instance variable for supporting the nextRecord() method
my $sorted = 0;

=item new()

B<Description:> Instantiate GFF3Builder object

B<Parameters:> 

 %args

B<Returns:> Returns a reference to GFF3Builder

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

    ## This lookup is keyed on GFF3 header name with value being the GFF3 header value.
    $self->{'_headers'} = {};

    ## This lookup is keyed on the id with value being a reference to the GFF3Record.
    ## Reference for every GFF3Record created by the GFF3Builder will be stored in
    ## this lookup.
    $self->{'_records'} = {};

    ## This lookup is keyed on an id (reference sequence) with value being a reference to an list of id values.
    ## The id values in the list are all of the features that are localized to the reference/master sequence.
    $self->{'_features'} = {};

    ## This lookup is keyed on the id with value being a reference to the GFF3Record.
    ## Unlike the _records lookup, this will only contain references to GFF3Records that
    ## have some FASTA sequence.
    $self->{'_recordswithfasta'} = {};

}

=item DESTROY

B<Description:> GFF3Builder class destructor

B<Parameters:> None

B<Returns:> None

=cut

sub DESTROY  {
    my $self = shift;

}

=item $obj->doesRecordExist($id)

B<Description:> Check whether builder has reference to some GFF3Record with $id

B<Parameters:> 

 $id - scalar/string corresponding to the ID attribute of the GFF3Record

B<Returns:> 

 0 - scalar false
 1 - scalar true

=cut

sub doesRecordExist {
    my ($self) = shift;
    my ($id) = @_;

    if (!defined($id)){
	$logger->logdie("id was not defined");
    }

    if (exists $self->{'_records'}->{$id}){
	return 1;
    }
    return 0;
}

=item $obj->addHeader($key, $value)

B<Description:> Add GFF3 file header

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

    if (! exists $self->{'_headers'}->{$key}){
      $self->{'_headers'}->{$key} = $value;
      push(@{$self->{'_headersArray'}}, "$key $value");
    } else {
      $logger->warn("header key '$key' value '$value' was already stored");
    }
}
 

=item $obj->getRecordById($id)

B<Description:> Return reference to GFF3Record with id $id

B<Parameters:> 

 $id - scalar/string value that corresponds to the ID attribute of the GFF3Record

B<Returns:> reference to the GFF3Record with id $id

=cut

sub getRecordById {

    my ($self) = shift;
    my ($id) = @_;

    if (!defined($id)){
	$logger->logdie("id was not defined");
    }
    
    if ( exists $self->{'_records'}->{$id} ){
	return $self->{'_records'}->{$id};
    }
    else {
	$logger->warn("GFF3Record does not exist for id '$id'");
    }

    return undef;
}
 
=item $obj->createAndAddRecord($id, $seqid, $source, $type, $start, $stop, $score, $strand, $phase)

B<Description:> Creates a new GFF3Record object and stores reference to that object

B<Parameters:> 

 $id     - the column 9 ID attribute
 $seqid  - column 1 seqid
 $source - column 2 source
 $type   - column 3 type
 $start  - column 4 start
 $stop   - column 5 stop
 $score  - column 6 score
 $strand - column 7 strand
 $phase  - column 8 phase

B<Returns:> reference to the GFF3Record

=cut

sub createAndAddRecord {

    my ($self) = shift;
    my ($id, $seqid, $source, $type, $start, $stop, $score, $strand, $phase) = @_;

    my $gff3Record = new GFF3::GFF3Record($id, $seqid, $source, $type, $start, $stop, $score, $strand, $phase, $self->{'_gbrowse'});
    if (!defined($gff3Record)){
	$logger->logdie("Could not instantiate GFF3::GFF3Record for id '$id' type '$type' start '$start' stop '$stop' strand '$strand'");
    }

    ## Store reference to this new GFF3Record
    $self->{'_records'}->{$id} = $gff3Record;

    ## Increment the class member
    $gff3RecordCtr++;

    ## Return reference to the GFF3Record
    return $gff3Record;
}

=item $obj->removeRecord($id)

B<Description:> Remove all references to the GFF3Record with id $id

B<Parameters:> 

 $id     - the column 9 ID attribute

B<Returns:> None

=cut

sub removeRecord {

    my ($self) = shift;
    my ($id) = @_;

    delete $self->{'_records'}->{$id};
    delete $self->{'_recordswithfasta'}->{$id};

    ## Keep track of which id values for GFF3Records that have been removed.
    $self->{'_removedrecords'}->{$id}++;
    ## Decrement the class member
    $gff3RecordCtr--;

}

=item $obj->extractFastaSequenceFromRecord($id)

B<Description:> Removes that FASTA sequence from the GFF3Record with id $id

B<Parameters:> 

 $id - scalar/string the ID attribute of the GFF3Record

B<Returns:> scalar/string FASTA sequence

=cut

sub extractFastaSequenceFromRecord {

    my ($self) = shift;
    my ($id) = @_;


    if ($self->doesRecordExist($id)){
	#    my $gff3Record = $self->{'_records'}->{$id};

	if ( exists $self->{'_recordswithfasta'}->{$id} ) {
	    ## Delete reference from the internal _recordswithfasta lookup
	    
	    delete $self->{'_recordswithfasta'}->{$id};
	}
	else {
	    $logger->warn("Trying to extract FASTA from GFF3Record with id '$id' ".
			  "but the GFF3Builder does not have reference to that record ".
			  "in the _recordswithfasta lookup");
	}

	return $self->{'_records'}->{$id}->extractFastaSequence();
    }
    else {
	$logger->logdie("There does not exist a GFF3Record with id '$id'");
    }

    return undef;
}

=item $obj->addFastaSequenceToRecord($id)

B<Description:> Adds FASTA sequence to the GFF3Record with id $id

B<Parameters:> 

 $id    - scalar/string the ID attribute of the GFF3Record
 $fasta - scalar/string the FASTA sequence

B<Returns:> None

=cut

sub addFastaSequenceToRecord {

    my ($self) = shift;
    my ($id, $fasta) = @_;

    if (!defined($id)){
	$logger->logdie("id was not defined");
    }
    if (!defined($fasta)){
	$logger->logdie("fasta was not defined for id '$id'");
    }

    if ($self->doesRecordExist($id)){
	## Retrieve reference to this GFF3Record
	my $gff3Record = $self->{'_records'}->{$id};

	## Add the FASTA sequence to that GFF3Record
	$gff3Record->addFastaSequence($fasta);

	## Add reference to the internal _recordswithfasta lookup
	$self->{'_recordswithfasta'}->{$id} = $gff3Record;
    }
    else {
	$logger->logdie("Attempting to add FASTA sequence to a GFF3Record that does not exist for id '$id'!");
    }
}

=item $obj->linkFeatureToSequence($sequence_id, $gff3Record)

B<Description:> Associates the GFF3Record of a feature to the contig GFF3Record

B<Parameters:> 

 $sequence_id - scalar/string the ID attribute of the contig GFF3Record
 $gff3Record  - GFF3Record

B<Returns:> None

=cut

sub linkFeatureToSequence {

    my ($self) = shift;
    my ($sequence_id, $gff3Record) = @_;

    if (!defined($sequence_id)){
	$logger->logdie("sequence_id was not defined");
    }

    if (!defined($gff3Record)){
	$logger->logdie("gff3Record was not defined");
    }

    my $feature_id = $gff3Record->getId();
    if (!defined($feature_id)){
	$logger->logdie("feature_id was not defined while.  sequence_id was '$sequence_id'");
    }
    
    push( @{$self->{'_features'}->{$sequence_id}}, $feature_id);
}

=item $obj->writeFile()

B<Description:> Calls methods that write the data to the GFF3 file defined by _filename

B<Parameters:> None

B<Returns:> None

=cut

sub writeFile {

    my ($self) = shift;

    if (exists $self->{'_filename'} ){
	my $fh;
	open ($fh, ">$self->{'_filename'}") || $logger->logdie("Could not open file '$self->{'_filename'}' in write mode: $!");
	$self->writeHeaders($fh);
	$self->writeRecords($fh);
	$self->writeFastaRecords($fh);
	close $fh;	
    }
    else {
	$logger->logdie("filename not defined");
    }
}

=item $obj->writeHeaders()

B<Description:> Writes the headers to the GFF3 file

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
	foreach my $header ( @{$self->{'_headersArray'}} ) {
	    print $fh "\#\#$header\n";
	}
    }
}

=item $obj->writeRecords()

B<Description:> Writes the GFF3Records to the GFF3 file

B<Parameters:> 

$fh - file handle

B<Returns:> None

=cut

sub writeRecords {

    my ($self) = shift;
    my ($fh) = @_;

    if (!defined($fh)){
	$logger->logdie("file handle was not defined");
    }
    
    foreach my $sequence_id ( keys %{$self->{'_features'}} ) {

	if (exists $self->{'_records'}->{$sequence_id}){
	    ## Write all of the primary sequences first.
	    $self->{'_records'}->{$sequence_id}->writeRecord($fh, $sequence_id, $self->{'_source'});
	    delete $self->{'_records'}->{$sequence_id};
	}
	else {
	    $logger->logdie("GFF3Record does not exist for id '$sequence_id'");
	}

	$self->_sortFeaturesByCoordinates($sequence_id);

	foreach my $feature_id ( @{$self->{'_features'}->{$sequence_id}} ) {
	    ## Write all of the features localized to the primary sequences next.
	    if (exists $self->{'_records'}->{$feature_id}){
		$self->{'_records'}->{$feature_id}->writeRecord($fh, $sequence_id, $self->{'_source'});
		delete $self->{'_records'}->{$feature_id};
	    }
	    else {
		if (exists $self->{'_removedrecords'}->{$feature_id} ){
		    if ($logger->is_debug()){
			## It could be that the GFF3Record was legitimately removed.
			## See the removeRecord() method.
			$logger->debug("GFF3Record with id '$feature_id' was removed.");
		    }
		}
		else {
		    $logger->logdie("GFF3Record with id '$feature_id' does not exist!");
		}			
	    }
	}
    }

    ## Write all of the primary sequences that did not have any features that
    ## were localized to the primary sequences.
    foreach my $id (keys %{$self->{'_records'}} ){
	$self->{'_records'}->{$id}->writeRecord($fh, $id, $self->{'_source'});
    }
}

=item $obj->writeFasta()

B<Description:> Writes the FASTA headers and sequences to the GFF3 file

B<Parameters:> 

$fh - file handle

B<Returns:> None

=cut

sub writeFastaRecords {
    my ($self) = shift;
    my ($fh) = @_;

    if (!defined($fh)){
	$logger->logdie("file handle was not defined");
    }

    if (exists $self->{'_recordswithfasta'}){

	print $fh "\#\#FASTA\n";

	foreach my $id ( keys %{$self->{'_recordswithfasta'}} ) {
	    $self->{'_recordswithfasta'}->{$id}->writeFasta($fh);
	}
    }
    else {
	$logger->warn("There weren't any GFF3Records that had FASTA sequence.");
    }
}

=item $obj->nextRecord()

B<Description:> Iteratively returns reference to each GFFRecord

B<Parameters:> None

B<Returns:> reference to GFF3Record

=cut

sub nextRecord {

    my ($self) = shift;

    if ($sorted == 0 ) {
	if ($logger->is_debug()){
	    $logger->debug("The GFF3Records have not yet been sorted by id");
	}
	## If the GFF3Records aren't already sorted, do so now.
	foreach my $id (sort keys %{$self->{'_records'}} ) {
	    push ( @{$self->{'_sortedrecords'}}, $id );	
	}

	## The list has been sorted once, let's not do it again.
	$sorted = 1;
    }

    if ( $recordIndex < $gff3RecordCtr ){
	my $id = $self->{'_sortedrecords'}->[$recordIndex++];
	if (defined($id)){
	    return $self->getRecordById($id);
	}
	else {
	    $logger->logdie("id was not defined for recordIndex '$recordIndex'");
	}
    }

    return undef;
}

=item $obj->_sortFeaturesByCoordinates()

B<Description:> Sort all of the GFFRecords by start coordinate

B<Parameters:> None

B<Returns:> None

=cut

sub _sortFeaturesByCoordinates {

    my $self = shift;
    my ($sequence_id) = @_;

    my $featuresList = $self->{'_features'}->{$sequence_id};
    my $recordsLookup = $self->{'_records'};

    my $sortedFeatureList=[];
    my $sortableFeatureList=[];
    use Data::Dumper;
    foreach my $feature_id ( @{$featuresList} ) {
	if (! exists $self->{'_records'}->{$feature_id}){
	    $logger->warn("feature_id '$feature_id' does not exist in the _records lookup!");
	}
	else {
	    push(@{$sortableFeatureList}, $feature_id);
	}
    }

    foreach my $feature_id ( sort { $self->{'_records'}->{$a}->getStart() <=> $self->{'_records'}->{$b}->getStart() }  @{$sortableFeatureList} ) {
	push(@{$sortedFeatureList}, $feature_id);
    }

    $self->{'_features'}->{$sequence_id} = $sortedFeatureList;
}

1; ## End of module
