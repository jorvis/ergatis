package GTF::GTFRecord;

=head1 NAME

GTF::GTFRecord.pm

=head1 VERSION

1.0

=head1 SYNOPSIS

use GTFRecord;

my $obj = new GTF::GTFRecord($id, $seqid, $source, $type, $start, $stop, $score, $strand, $phase);

$obj->addParent($parentId);

$obj->addFastaSequence("GTCGTAGCTCAGCATGA");

$obj->writeRecord($fileHandle1);

=head1 AUTHOR

Jay Sundaram
sundaram@tigr.org

=head1 METHODS


=over 4

=cut

use strict;
use GTF::GTFBuilder;
use GTF::Logger;

my $logger = GTF::Logger::get_logger("Logger::GTF");


## Class variable for keeping track of the number of GTF records
my $gtfRecordCtr = 0;

=item new($seqname, $source, $feature, $start, $end, $score, $frame, $attributes, $comments)

B<Description:> Instantiate GTFRecord object

B<Parameters:> 

 $seqname     - column 1
 $source      - column 2
 $feature     - column 3
 $start       - column 4
 $end         - column 5
 $score       - column 6
 $strand      - column 7
 $frame       - column 8
 $attributes  - column 9
 $comments    - column 10

B<Returns:> reference to the GTFRecord

=cut

sub new  {

    my $class = shift;


    my ($seqname, $source, $feature, $start, $end, $score, $strand,
	$frame, $attributes, $comment) = @_;


   my $self = {};


    if (!defined($seqname)){
	$logger->logdie("seqname was not defined");
    }
    else {
	$self->{'_core'}->[0] = $seqname;
	$self->{'_seqname'} = $seqname;
    }

    if (defined($source)){
	$self->{'_core'}->[1] = $source;
    }

    if (defined($feature)){
	$self->{'_core'}->[2] = $feature;
    }
    
    if (!defined($start)){
	$logger->logdie("start was not defined");
    }
    else {
	$self->{'_core'}->[3] = $start;
    }

    if (!defined($end)){
	$logger->logdie("end was not defined");
    }
    else {
	$self->{'_core'}->[4] = $end;
    }

    if (!defined($score)){
	$logger->logdie("score was not defined");
    }
    else {
	$self->{'_core'}->[5] = $score;
    }

    if (!defined($strand)){
	$strand = '+';
    }
    $self->{'_core'}->[6] = $strand;

    if (!defined($frame)){
	$frame = 0;
    }
    $self->{'_core'}->[7] = $frame;



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
}

=item DESTROY

B<Description:> GTFRecord class destructor

B<Parameters:> None

B<Returns:> None

=cut 

sub DESTROY  {

    my $self = shift;

}

=item $obj->addAttribute($key, $value)

B<Description:> Add named attributes to the GTFRecord object

B<Parameters:>

 $key   - name of the attribute
 $value - value to add to the record

B<Returns:>  None

=cut 

sub addAttribute {

    my $self = shift;
    my ($key, $value) = @_;

    if (!defined($key)){
	$logger->logdie("key was not defined");
    }
    if (!defined($value)){
	$logger->logdie("value was not defined");
    }

    ## URL encoding for compliance with GTF spec
    $value =~ s/\\n/\n/g;
    $value =~ s/([\f\n\r\t;=%&,[:cntrl:]])/sprintf("%%%02X", ord($1))/seg;

    push(@{$self->{'_attrs'}->{$key}}, $value);

    if ( ($key eq 'gene_product_name') || ($key eq 'description')){
      if ((exists $self->{'_gbrowse'}) && (defined($self->{'_gbrowse'})) && ($self->{'_gbrowse'} == 1)){
	## In order to ensure GBrowse compatibility, the gene_product_name
	## must be store as the Note attribute
	push(@{$self->{'_attrs'}->{'Note'}}, $value);
      }
    }
}

=item $obj->addParent($parent)

B<Description:> Add Parent attribute to the GTFRecord object

B<Parameters:>

$parent - string identifier that matches the ID attribute for the parent GTFRecord

B<Returns:> None

=cut 

sub addParent {

    my $self = shift;
    my ($parent) = @_;

    if (!defined($parent)){
	$logger->logdie("parent was not defined");
    }

    $self->{'_parent'} = $parent;
}

=item $obj->getParent()

B<Description:> Retrieves the Parent attribute for the GTFRecord object

B<Parameters:> None

B<Returns:>  Scalar/string that is the value assigned to the Parent attribute for the GTFRecord

=cut 

sub getParent {

    my $self = shift;

    if ( exists $self->{'_parent'} ) {
	return $self->{'_parent'};
    }
    
    return undef;
}

=item $obj->getId()

B<Description:> Retrieves the ID attribute for the GTFRecord object

B<Parameters:> None

B<Returns:> Scalar/string that is the ID assigned to the GTFRecord

=cut 

sub getId {

    my $self = shift;
    
    if ( exists $self->{'_id'} ) {

	return  $self->{'_id'};

    }  else {
	$logger->logdie("GTFRecord with no id!");
    }
}

=item $obj->getType()

B<Description:> Retrieves the type for the GTFRecord object

B<Parameters:> None

B<Returns:> Scalar/string that is the type assigned to the GTFRecord

=cut 

sub getType {

    my $self = shift;

    if ( exists $self->{'_core'} ) {

	if (defined($self->{'_core'}->[2])) {

	    return  $self->{'_core'}->[2];
	}
	else {
	    $logger->logdie("The type is not defined for the ".
			    "GTFRecord with id '$self->{'_id'}'");
	}

    } else {
	$logger->logdie("The _core data does not exist for the ".
			"GTFRecord with id '$self->{'_id'}'!");
    }
}

=item $obj->getStart()

B<Description:> Retrieves the start for the GTFRecord object

B<Parameters:> None

B<Returns:> Scalar/string that is the start assigned to the GTFRecord

=cut 

sub getStart {

    my $self = shift;

    if ( exists $self->{'_core'} ) {

	if (defined($self->{'_core'}->[3])) {

	    return  $self->{'_core'}->[3];

	} else {

	    $logger->logdie("The start is not defined for the ".
			    "GTFRecord with id '$self->{'_id'}'");
	}

    } else {

	$logger->logdie("The _core data does not exist for the ".
			"GTFRecord with id '$self->{'_id'}'!");
    }
}

=item $obj->extractFastaSequence()

B<Description:> Returns the FASTA sequence (scalar/string) and deletes the fasta attribute for the GTFRecord object

B<Parameters:> None

B<Returns:> Scalar/string that is the FASTA sequence belonging to the GTFRecord

=cut 
 
sub extractFastaSequence {

    my $self = shift;

    if ( exists $self->{'fasta'} ){

	my $fasta = $self->{'fasta'};

	delete $self->{'fasta'};

	return $fasta;

    } else {

	$logger->warn("No FASTA for GTFRecord with ".
		      "id '$self->{'_id'}'");
    }

    return undef;
}

=item $obj->hasFasta()

B<Description:> Checks whether a FASTA sequence has been assigned to the GTFRecord object

B<Parameters:> None

B<Returns:> 

 1 - true
 0 - false

=cut 

sub hasFasta {

    my $self = shift;

    if ( exists $self->{'fasta'} ){
	return 1;
    }

    return 0;
}

=item $obj->addFastaSequence($fasta)

B<Description:>  Adds a FASTA sequence to the GTFRecord object

B<Parameters:>  Scalar/string that is the FASTA sequence

B<Returns:> None

=cut 

sub addFastaSequence {

    my $self = shift;
    my ($fasta) = @_;

    if (!defined($fasta)){
	$logger->logdie("fasta was not defined");
    }

    $self->{'fasta'} = $fasta;
}
 
=item $obj->writeRecord($fh, $seqid, $source)

B<Description:> Method that writes the GTF record to the filehandle fh.

B<Parameters:>

 $fh     - filehandle for output file
 $seqid  - column 1 seqid
 $source - column 2 source

B<Returns:>  None

=cut 

sub writeRecord {

    my $self = shift;
    my ($fh, $seqid, $source) = @_;


    if (!defined($fh)){
	$logger->logdie("fh was not defined");
    }
    if (!defined($seqid)){
	$logger->logdie("seqid was not defined");
    }
    if (!defined($source)){
	$logger->logdie("source was not defined");
    }

    my $cols = $self->{'_core'};

    if (defined($cols->[0])){
	$seqid = $cols->[0];
    }

    if (defined($cols->[1])){
	$source = $cols->[1];
    }

    ##         seqid   source   type        start       stop        score       strand      phase      
    print $fh "$seqid\t$source\t$cols->[2]\t$cols->[3]\t$cols->[4]\t$cols->[5]\t$cols->[6]\t$cols->[7]\t";

    my $values = "ID=$self->{'_id'};";

    my $parent = $self->getParent();
    if (defined($parent)){
	$values .= "Parent=$parent;";
    }

    my $nameAttributeEncountered=0;

    foreach my $attr ( keys %{$self->{'_attrs'}} ) {

      if ($attr eq 'Name'){
	$nameAttributeEncountered++;
      }

	$values .= "$attr=";

	foreach my $value ( @{$self->{'_attrs'}->{$attr}} ) {
	    $values .= "$value,";
	}
	## remove trailing comma
	chop $values;
	
	## add attribute separator (semicolon)
	$values .= ";";
    }

    if ( ($nameAttributeEncountered == 0) && (($cols->[2] eq 'mRNA') || ( $cols->[2] eq 'contig' ) || ( $cols->[2] eq 'assembly' )) ){
      ## Need to add the Name attribute for records
      ## of type mRNA, contig or assembly
      ## in order for GTF file to be valid.
      $values .= "Name=$self->{'_id'};";
    }


    ## remove trailing semicolon
    chop $values;

    print $fh "$values\n";
}

=item $obj->writeFasta($fh)

B<Description:> Method that writes the FASTA sequence for the GTFRecord to the filehandle fh.

B<Parameters:>

 $fh - filehandle for output file

B<Returns:> None

=cut 

sub writeFasta {

    my $self = shift;
    my ($fh) = @_;

    if (!defined($fh)){
	$logger->logdie("fh was not defined");
    }

    if ( $self->hasFasta() ){
	my $formattedFasta = &formatFasta($self->{'_id'}, $self->{'fasta'});
	print $fh "$formattedFasta";
    }
    else {
	$logger->fatal("GTFRecord with id '$self->{'_id'}' does not have FASTA");
    }
}

=item $obj->formatFasta($fastaHeader, $fastaSequence)

B<Description:> Formats the FASTA header and sequence

B<Parameters:>

 $fastaHeader   - scalar/string
 $fastaSequence - scalar/string

B<Returns:> $fastaRecord - scalar/string formatted GTF FASTA record

=cut 

sub formatFasta {

    my ($fastaHeader, $fastaSequence) = @_;
    #This subroutine takes a sequence name and its sequence and
    #outputs a correctly formatted single fasta entry (including newlines).
    
    my $fastaRecord=">"."$fastaHeader"."\n";

    $fastaSequence =~ s/\s+//g; ## remove white spaces

    for(my $i=0; $i < length($fastaSequence); $i+=60){

	my $seq_fragment = substr($fastaSequence, $i, 60);

	$fastaRecord .= "$seq_fragment"."\n";
    }
    return $fastaRecord;

}

=item $obj->getTranslationTable()

B<Description:> Returns the translation_table value

B<Parameters:> None

B<Returns:>  scalar/string translation_table value

=cut 

sub getTranslationTable {

    my $self = shift;

    if (exists $self->{'_attrs'}->{'translation_table'}){

	if (defined($self->{'_attrs'}->{'translation_table'}->[0])){

	    return $self->{'_attrs'}->{'translation_table'}->[0];

	} else {

	    $logger->warn("translation_table is not defined");
	}

    } else {

	$logger->warn("translation_table attribute does not exist");
    }

    return undef;
}

=item $obj->addTranslationTable()

B<Description:> Stores a translation_table attribute for the GTFRecord

B<Parameters:> 

 $translation_table - scalar/string value that should be an integer value representing the translation_table

B<Returns:>  None

=cut 

sub addTranslationTable {

    my $self = shift;
    my ($translation_table) = @_;

    if (!defined($translation_table)){
	$logger->logdie("translation_table was not defined");
    }

    push(@{$self->{'_attrs'}->{'translation_table'}}, $translation_table);
}

=item $obj->hasTranslationTable()

B<Description:> Verifies whether the translation_table attribute has been defined

B<Parameters:> None

B<Returns:>

 0 - false
 1 - true

=cut 

sub hasTranslationTable {
    
    my $self = shift;
    
    if (exists $self->{'_attrs'}->{'translation_table'}){

	if (defined($self->{'_attrs'}->{'translation_table'}->[0])){

	    return 1;
	}
    }

    return 0;
}

1==1; ## end of module
