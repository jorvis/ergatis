package GTF::GTFBuilder;

=head1 NAME

GTF::GTFBuilder.pm 

A class to facilitate the creation of GTF files.

Will be based on the specification outlined here:
http://mblab.wustl.edu/GTF22.html#resources


=head1 VERSION

1.0

=head1 SYNOPSIS

use GTFBuilder;

my $obj = new GTF::GTFBuilder( filename => '/tmp/myfile.gtf');

$obj->addRecord($gtfRecord);

$obj->writeRecords();

=head1 AUTHOR

Jay Sundaram
sundaram@tigr.org

=head1 METHODS


=over 4

=cut

use strict;
use GTF::GTFRecord;
use GTF::Logger;
use Carp;

my $logger = GTF::Logger::get_logger("Logger::GTF");


## Class variable for keeping track of the number of GTF records
my $gtfRecordCtr = 0;

## Class variable for keeping track of the number of GTF records with FASTA sequence.
my $fastaRecordCtr = 0;

## Instance variable for supporting the nextRecord() method
my $recordIndex = 0;

## Instance variable for supporting the nextRecord() method
my $sorted = 0;

=item new()

B<Description:> Instantiate GTFBuilder object

B<Parameters:> 

 %args

B<Returns:> Returns a reference to GTFBuilder

=cut

sub new  {

    my $class = shift;
    my $self = {};
    bless $self, $class;    

    $self->_init(@_);

    confess "Not fully implmented/tested";

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

    ## This lookup is keyed on GTF header name with value being the GTF header value.
    $self->{'_headers'} = {};

    ## This lookup is keyed on the id with value being a reference to the GTFRecord.
    ## Reference for every GTFRecord created by the GTFBuilder will be stored in
    ## this lookup.
    $self->{'_records'} = {};

    ## This lookup is keyed on an id (reference sequence) with value being a reference to an list of id values.
    ## The id values in the list are all of the features that are localized to the reference/master sequence.
    $self->{'_features'} = {};

    ## This lookup is keyed on the id with value being a reference to the GTFRecord.
    ## Unlike the _records lookup, this will only contain references to GTFRecords that
    ## have some FASTA sequence.
    $self->{'_recordswithfasta'} = {};

}

=item DESTROY

B<Description:> GTFBuilder class destructor

B<Parameters:> None

B<Returns:> None

=cut

sub DESTROY  {
    my $self = shift;

}

=item $obj->doesRecordExist($id)

B<Description:> Check whether builder has reference to some GTFRecord with $id

B<Parameters:> 

 $id - scalar/string corresponding to the ID attribute of the GTFRecord

B<Returns:> 

 0 - scalar false
 1 - scalar true

=cut

sub doesRecordExist {

    my $self = shift;
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

B<Description:> Add GTF file header

B<Parameters:> 

 $key   - scalar/string name of the header
 $value - scalar/string value of the header attribute

B<Returns:> None

=cut

sub addHeader {

    my $self = shift;
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

B<Description:> Return reference to GTFRecord with id $id

B<Parameters:> 

 $id - scalar/string value that corresponds to the ID attribute of the GTFRecord

B<Returns:> reference to the GTFRecord with id $id

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
	$logger->warn("GTFRecord does not exist for id '$id'");
    }

    return undef;
}
 
=item $obj->createAndAddRecord($id, $seqid, $source, $type, $start, $stop, $score, $strand, $phase)

B<Description:> Creates a new GTFRecord object and stores reference to that object

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

B<Returns:> reference to the GTFRecord

=cut

sub createAndAddRecord {

    my $self = shift;
    my ($id, $seqid, $source, $type, $start, $stop, $score, $strand, $phase) = @_;

    my $gtfRecord = new GTF::GTFRecord($id, $seqid, $source, $type, $start, $stop, $score, $strand, $phase, $self->{'_gbrowse'});
    if (!defined($gtfRecord)){
	$logger->logdie("Could not instantiate GTF::GTFRecord for id '$id' type '$type' start '$start' stop '$stop' strand '$strand'");
    }

    ## Store reference to this new GTFRecord
    $self->{'_records'}->{$id} = $gtfRecord;

    ## Increment the class member
    $gtfRecordCtr++;

    ## Return reference to the GTFRecord
    return $gtfRecord;
}

=item $obj->removeRecord($id)

B<Description:> Remove all references to the GTFRecord with id $id

B<Parameters:> 

 $id     - the column 9 ID attribute

B<Returns:> None

=cut

sub removeRecord {

    my $self = shift;
    my ($id) = @_;

    delete $self->{'_records'}->{$id};
    delete $self->{'_recordswithfasta'}->{$id};

    ## Keep track of which id values for GTFRecords that have been removed.
    $self->{'_removedrecords'}->{$id}++;
    ## Decrement the class member
    $gtfRecordCtr--;

}

=item $obj->extractFastaSequenceFromRecord($id)

B<Description:> Removes that FASTA sequence from the GTFRecord with id $id

B<Parameters:> 

 $id - scalar/string the ID attribute of the GTFRecord

B<Returns:> scalar/string FASTA sequence

=cut

sub extractFastaSequenceFromRecord {

    my $self = shift;
    my ($id) = @_;


    if ($self->doesRecordExist($id)){
	#    my $gtfRecord = $self->{'_records'}->{$id};

	if ( exists $self->{'_recordswithfasta'}->{$id} ) {
	    ## Delete reference from the internal _recordswithfasta lookup
	    
	    delete $self->{'_recordswithfasta'}->{$id};
	}
	else {
	    $logger->warn("Trying to extract FASTA from GTFRecord with id '$id' ".
			  "but the GTFBuilder does not have reference to that record ".
			  "in the _recordswithfasta lookup");
	}

	return $self->{'_records'}->{$id}->extractFastaSequence();
    }
    else {
	$logger->logdie("There does not exist a GTFRecord with id '$id'");
    }

    return undef;
}

=item $obj->addFastaSequenceToRecord($id)

B<Description:> Adds FASTA sequence to the GTFRecord with id $id

B<Parameters:> 

 $id    - scalar/string the ID attribute of the GTFRecord
 $fasta - scalar/string the FASTA sequence

B<Returns:> None

=cut

sub addFastaSequenceToRecord {

    my $self = shift;
    my ($id, $fasta) = @_;

    if (!defined($id)){
	$logger->logdie("id was not defined");
    }
    if (!defined($fasta)){
	$logger->logdie("fasta was not defined for id '$id'");
    }

    if ($self->doesRecordExist($id)){
	## Retrieve reference to this GTFRecord
	my $gtfRecord = $self->{'_records'}->{$id};

	## Add the FASTA sequence to that GTFRecord
	$gtfRecord->addFastaSequence($fasta);

	## Add reference to the internal _recordswithfasta lookup
	$self->{'_recordswithfasta'}->{$id} = $gtfRecord;
    }
    else {
	$logger->logdie("Attempting to add FASTA sequence to a GTFRecord that does not exist for id '$id'!");
    }
}

=item $obj->linkFeatureToSequence($sequence_id, $gtfRecord)

B<Description:> Associates the GTFRecord of a feature to the contig GTFRecord

B<Parameters:> 

 $sequence_id - scalar/string the ID attribute of the contig GTFRecord
 $gtfRecord  - GTFRecord

B<Returns:> None

=cut

sub linkFeatureToSequence {

    my $self = shift;
    my ($sequence_id, $gtfRecord) = @_;

    if (!defined($sequence_id)){
	$logger->logdie("sequence_id was not defined");
    }

    if (!defined($gtfRecord)){
	$logger->logdie("gtfRecord was not defined");
    }

    my $feature_id = $gtfRecord->getId();
    if (!defined($feature_id)){
	$logger->logdie("feature_id was not defined while.  sequence_id was '$sequence_id'");
    }
    
    push( @{$self->{'_features'}->{$sequence_id}}, $feature_id);
}

=item $obj->writeFile()

B<Description:> Calls methods that write the data to the GTF file defined by _filename

B<Parameters:> None

B<Returns:> None

=cut

sub writeFile {

    my $self = shift;

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

B<Description:> Writes the headers to the GTF file

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

B<Description:> Writes the GTFRecords to the GTF file

B<Parameters:> 

$fh - file handle

B<Returns:> None

=cut

sub writeRecords {

    my $self = shift;
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
	    $logger->logdie("GTFRecord does not exist for id '$sequence_id'");
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
			## It could be that the GTFRecord was legitimately removed.
			## See the removeRecord() method.
			$logger->debug("GTFRecord with id '$feature_id' was removed.");
		    }
		}
		else {
		    $logger->logdie("GTFRecord with id '$feature_id' does not exist!");
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

B<Description:> Writes the FASTA headers and sequences to the GTF file

B<Parameters:> 

$fh - file handle

B<Returns:> None

=cut

sub writeFastaRecords {
    my $self = shift;
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
	$logger->warn("There weren't any GTFRecords that had FASTA sequence.");
    }
}

=item $obj->nextRecord()

B<Description:> Iteratively returns reference to each GFFRecord

B<Parameters:> None

B<Returns:> reference to GTFRecord

=cut

sub nextRecord {

    my $self = shift;

    if ($sorted == 0 ) {
	if ($logger->is_debug()){
	    $logger->debug("The GTFRecords have not yet been sorted by id");
	}
	## If the GTFRecords aren't already sorted, do so now.
	foreach my $id (sort keys %{$self->{'_records'}} ) {
	    push ( @{$self->{'_sortedrecords'}}, $id );	
	}

	## The list has been sorted once, let's not do it again.
	$sorted = 1;
    }

    if ( $recordIndex < $gtfRecordCtr ){
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

1==1; ## end of module
