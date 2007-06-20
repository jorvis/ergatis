package Chado::Gene;

=head1 NAME

  Gene.pm - Object representation of a gene

=head1 SYNOPSIS

use Gene;

my $gene = new Gene('geneID');

#addWhatever(id, fmin, fmax, strand(0,1)), 
$gene->addExon('exonID', 5, 100, 0);
$gene->addPolypeptide('polypeptideID', 5, 100, 0);
$gene->addTranscript('transcriptID', 5, 100, 0);
$gene->addCDS('cdsID', 5, 100, 0);
$gene->addFeature('startCodonID', 0, 3, 0);

#Create a group (each group represents an alternate form)
#(ie, alternative splicing).
$filter = { 'all' => 1 };
$gene->createGroup('groupID', $filter); 

=head1 ADVANCED USAGE

    For adding features to groups, you can use a filter hash reference which
    allows the user to filter the types, id, or strand of features to include in
    a group.  The format would be as follows:

    $filter = { 'id' => 'exonID' };
    $filter = { 'type' => 'exon' };
    
    Or if you wanted all features in the gene to be added to a certain group
    you could use the all key with any non-zero value.

    $filter = { 'all' => 1 };

    When filtering on more than one parameter, the functionality defaults to
    include one of the listed parameters.  For example:

    $filter = { 'type' => 'exon',
                'strand' => 0 };

    Using this filter would include all exons and all features on the forward
    strand.  If instead you wanted only exons on the forward strand, use the 
    'and' key:

    $filter = { 'type' => 'exon',
                'strand' => 0,
                'and' => 1 };

    You can also filter on the start and stop coordinates.  To indicated greater than a plus
    should precede the number (be sure to put value in quotes to preserve the plus/minus sign).
    To indicate less than, you guessed it, a minus sign.  Example:

    $filter = { 'start' => "+6000" };

    Would include those features that had a start position of 6000 or higher.  

    The filter can be used in the createGroup, addToGroup, or accessor (get) 
    subroutines.  

    

    NOTE:
    Currently, there are no checks to ensure the inclusion of certain features
    within a gene.  For example, a gene doesn't have to have a CDS, polypeptide,
    or transcript feature to be a valid gene.

=head1 CONTACT

    Kevin Galens
    kgalens@tigr.org

=cut

use strict;
use warnings;
use Class::Struct;

#######A Data Structure Definition####################
 struct( Feature => {
        'id'        => '$',
        'start'     => '$',
        'stop'      => '$',
        'strand'    => '$',
        'type'      => '$',
        'attribute' => '*%',
        'children'  => '*%',  #<-- Example of how we can store relationships.
    });
##########################################################


############Constructor############################################
=item new Gene( $id, $start, $stop $strand, $seqId, $score );

B<Description:> Creates a gene object with an id.

B<Parameters:>  
    $id - a unique identifier for the gene
    $start - start nucleotide of the gene
    $stop - end nucleotide of the gene
    $strand - 0 (forward) or 1 (reverse)
    $seqId - an identifier for the nucleotide sequence on which the gene is located 
    $score - the score of the gene (if predicted).

B<Returns:> a gene object reference

=cut 

sub new {
    my $class = shift;
    my $self = {};

    #Bless it as a gene.
    bless $self, $class;

    #Initialize the data type.
    $self->_init(@_);

    #And Return
    return $self;
}

###########################CDS Manipulation################################

=item $gene->addCDS( $id, $fmin, $fmax, $strand );

B<Description:> Adds one CDS to the gene.

B<Parameters:> $id = unique id of the CDS
               $fmin = start boundary of the CDS
               $fmax = end boundary of the CDS
               $strand = forward (0) or reverse (1)

B<Returns:> The id of the CDS added.

=cut

sub addCDS {
    my ($self, @properties) = @_;
    return $self->addFeature(@properties, 'cds');
}

=item $gene->addMultiCDS($CDSs);

B<Description:> Sets all the CDSs of the gene.  Will remove all previously associated
    CDSs from the gene.

B<Parameters:> Array reference of Gene::Feature objects.

B<Returns:> Number of CDSs added.

=cut

sub addMultiCDS {
    #Sets the coding sequence for the gene
    my ($self, $cds) = @_;
    return $self->addMultiFeatures($cds, 'cds');
}

=item $gene->getCDS([$filter]);

B<Description:> Returns the CDSs related to the gene object based on
    the filter.

B<Parameters:> Optional filter parameter (see Gene.pm perldoc for 
                                          format of hash ref)

B<Returns:> Reference to an array of cds.  Format:
    [$cds1, $cds2,...] 

    where 
    $cds->{(id|start|stop|strand|type}} = value

=cut

sub getCDS {
    my ($self, $filter) = @_;
    $filter = { 'all' => 1 } unless($filter);
    return $self->getFeatures($filter, 'cds');
}

#######################Exon Manipulation################################

=item $gene->addExon( $id, $fmin, $fmax, $strand );

B<Description:> Adds one exon to the gene.

B<Parameters:> $id = unique id of the exon
               $fmin = start boundary of the exon
               $fmax = end boundary of the exon
               $strand = forward (0) or reverse (1)

B<Returns:> The id of the exon added.

=cut

sub addExon {
    #Adds an exon
    my ($self, @properties) = @_;
    return $self->addFeature(@properties, 'exon');
}

=item $gene->addMultiExons($exons);

B<Description:> Sets all the exons of the gene.  Will remove all exons
    that were previously associated with gene.

B<Parameters:> Array Reference of Gene::Feature objects

B<Returns:> Number of exons added.

=cut

sub addMultiExons {
    #Sets the exons of the gene.  All of them
    my ($self, $exons) = @_;
    return $self->addMultiExons($exons, 'exon');
}

=item $gene->getExons([$filter]);

B<Description:> Returns the exons related to the gene object based on
    the filter.

B<Parameters:> Optional filter parameter (see Gene.pm perldoc for 
                                          format of hash ref)

B<Returns:> Reference to an array of cds.  Format:
    [$cds1, $cds2,...] 

    where 
    $cds->{(id|start|stop|strand|type}} = value

=cut

sub getExons {
    my ($self, $filter) = @_;
    $filter = {'all'=>1} unless($filter);
    return $self->getFeatures($filter, 'exon');
}

#################Feature manipulation###############################

=item $gene->addFeature( $id, $fmin, $fmax, $strand, $type );

B<Description:> Adds one feature to the gene.

B<Parameters:> $id = unique id of the exon
               $fmin = start boundary of the exon
               $fmax = end boundary of the exon
               $strand = forward (0) or reverse (1)
               $type = the type of feature added (should be a SO term).

B<Returns:> The id of the feature added.

=cut

sub addFeature {
    #Adds a general feature to a gene (ex Poly A, start_codon etc.)
    my ($self, $id, $start, $stop, $strand, $type) = @_;
    my $feat;

    if($id->isa('Gene::Feature')) {
        $feat = $id; 
        $id = $feat->id;
    } else {
        
        #Create a new feature
        $feat = Feature->new(
                             id     => $id,
                             start  => $start,
                             stop   => $stop,
                             strand => $strand,
                             type   => $type,
                            );
    }
                          
    push(@{$self->{'features'}}, $feat);

    return $id;

}

=item $gene->addMultiFeatures($features);

B<Description:> Sets multiple features of the gene.  

B<Parameters:> Array reference of features

B<Returns:> Number of features added.

=cut

sub addMultiFeatures {
    #Sets all the features
    my ($self, $features, $type) = @_;

    my $startNum = scalar @{$self->{'features'}};
    foreach(@{$features}) {
        my $id = $self->addFeatureObj($_, $type);
    }
    my $endNum = scalar @{$self->{'features'}};

    return ($endNum-$startNum);
    
}

=item $gene->getFeatures([$filter]);

B<Description:> Returns the features related to the gene object based on
    the filter.

B<Parameters:> Optional filter parameter (see Gene.pm perldoc for 
                                          format of hash ref)

B<Returns:> Returns an array reference of features.
    $feature->{(id|start|stop|strand|type)}

=cut

sub getFeatures {
    #Returns the features
    my ($self,$filter,$type) = @_;
    $filter = {'all' => 1} unless($filter);
    $filter = &_mergHashRef($filter, { 'type' => $type }) if($type);
    my $retval = $self->filterFunction($filter);
    return $retval;
}

=item $gene->getFeature($id);

B<Description:> Returns a single feature reference for the feature with specified id.

B<Parameters:> feature id

B<Returns:> Returns a feature reference    

=cut

sub getFeature {
    my ($self, $featId) = @_;

    my $filter = { 'id' => $featId };
    my $feature = $self->filterFunction( $filter );

    if (scalar(@{$feature}) > 0) {
        return $feature->[0];
    } else {
        return undef;
    }
}

=item $gene->addFeatureScore($featId, $scoreType, $value);

B<Description:> Sets the score type for a feature

B<Parameters:> $featId - The id of the feature to associate the score with
               $scoreType - The type of score (ex. p-value, bit_score, etc.)
                            Should be found in an ontology
               $value - The actual score.

B<Returns:> Returns the id of the feature to which the score was added.
    
=cut

sub addFeatureScore {
    my ($self, $featId, $scoreType, $value) = @_;
    my $filter = { 'id' => $featId };
    my $feature_array = $self->filterFunction( $filter );
    my $feature = shift( @{$feature_array} );
    $feature->attribute($scoreType, $value);
    return $feature->id();
}

=item $gene->addFeatureAttribute($featId, $name, $content);

B<Description:> Adds an attribute under a feature.

B<Parameters:> $featId - The id of the feature to add an attribute to
               $name - Attribute name
               $content - Attribute content

B<Returns:> Returns feature reference for the modified feature or undef if id was not found
    
=cut

sub addFeatureAttribute {
    my ($self, $featId, $name, $content) = @_;
    
    my $feature = $self->getFeature($featId);

    if (defined($feature)) {
        $feature->{'Feature::attribute'}->{$name} = $content;
        return $feature;
    } else {
        return undef;
    }
}

#################ID manipulation#############################################
=item $gene->getId;

B<Description:> Get the id of the gene object

B<Parameters:> None

B<Returns:> Scalar id of the gene object.

=cut

sub getId {
    my $self = shift;
    return $self->{'id'};
}

=item $gene->setId( $newID );

B<Description:> Set the id of the gene object

B<Parameters:> Scalar id of the gene object.

B<Returns:> The new ID of the gene.

=cut

sub setId {
    my ($self,$id) = @_;
    $self->{'id'} = $id;

    return $self->{'id'};
}


#############################Polypeptide Manipulation#######################
=item $gene->addPolypeptide( $id, $fmin, $fmax, $strand );

B<Description:> Adds one polypeptide to the gene.

B<Parameters:> $id = unique id of the polypeptide
               $fmin = start boundary of the polypeptide
               $fmax = end boundary of the polypeptide
               $strand = forward (0) or reverse (1)

B<Returns:> The id of the polypeptide added.

=cut

sub addPolypeptide {
    my ($self, @properties) = @_;
    return $self->addFeature(@properties, 'polypeptide');
}

=item $gene->addMulitPolypeptide($polypeptides);

B<Description:> Sets all the polypeptides of the gene.  Will remove all polypeptides
    that were previously associated with gene.

B<Parameters:> Array Reference of Gene::Feature objects

B<Returns:> Number of polypeptides added.

=cut

sub addMultiPolypeptide {
    #Set the polypeptide
    my ($self, $polypeptides) = @_;
    return $self->addMultiFeatures($polypeptides, 'polypeptide');
}

=item $gene->getPolypeptide([$filter]);

B<Description:> Returns the polypeptides related to the gene object based on
    the filter.

B<Parameters:> Optional filter parameter (see Gene.pm perldoc for 
                                          format of hash ref)

B<Returns:> Hash Reference.  Format:
    $hashRef->{polypeptideID}->{(start|stop|strand|type)} = value

=cut

sub getPolypeptide {
    #Returns the polypeptide
    my ($self, $filter) = @_;
    $filter = {'all' => 1} unless($filter);
    return $self->getFeatures($filter, 'polypeptide');
}

#######################Transcript Manipulation#####################

=item $gene->addTranscript( $id, $fmin, $fmax, $strand );

B<Description:> Adds one transcript to the gene.

B<Parameters:> $id = unique id of the transcript
               $fmin = start boundary of the transcript
               $fmax = end boundary of the transcript
               $strand = forward (0) or reverse (1)

B<Returns:> The id of the transcript added.

=cut

sub addTranscript {
    my ($self, @properties) = @_;
    return $self->addFeature(@properties, 'transcript');
}


=item $gene->addMultiTranscript($transcripts);

B<Description:> Sets all the transcripts of the gene.  Will remove all transcripts
    that were previously associated with gene.

B<Parameters:> Reference to an array of Gene::Feature objects

B<Returns:> Number of transcripts added.

=cut

sub addMultiTranscript {
    my ($self, $transcripts) = @_;
    return $self->addMultiFeatures($transcripts, 'transcript');
}

=item $gene->getTranscript([$filter]);

B<Description:> Returns the transcripts related to the gene object based on
    the filter.

B<Parameters:> Optional filter parameter (see Gene.pm perldoc for 
                                          format of hash ref)

B<Returns:> Hash Reference.  Format:
    $hashRef->{transcriptID}->{(start|stop|strand|type)} = value

=cut

sub getTranscript {
    #Gets the transcript
    my ($self, $filter) = @_;
    $filter = {'all'=>1};
    return $self->getFeatures($filter, 'transcript');
}

#########################Sequence Id Manipulation#########################
=item $gene->getSeqId;

B<Description:> Returns the sequence upon which the gene is located.

B<Parameters:> None

B<Returns:> The sequence id

=cut

sub getSeqId {
    my $self = shift;
    return $self->{'seq'};
}

=item $gene->setSeqId($seqId);

B<Description:> Sets sequence upon which the gene is located.

B<Parameters:> The new sequence id

B<Returns:> The new sequence id

=cut

sub setSeqId {
    my ($self, $seqId) = @_;
    $self->{'seq'} = $seqId;
    return $seqId;
}

#######################Group Functions ###################################
=item $gene->addToGroup($groupID[, $filter]);

B<Description:> Add a pre-existing feature to specified group identified by 
    group id.  If the group does not exist, it will be created.  If filter is
    not passed, all features will be added to the group.

B<Parameters:> $groupID - Identifier of the group
               $filter - see format in Gene.pm perldoc (hash ref)

B<Returns:> The number of features added to the specified group.

=cut

sub addToGroup {
    #Adds to specified group
    my ($self, $groupID, $filter) = @_;
    my $results = $self->{'features'};
    $results = $self->filterFunction($filter) if($filter);
    my $retval = scalar @{$results};

    foreach my $feat (@{$results}) {
        $self->{'groups'}->{$groupID}->{$feat->id} = $feat;
    }

    return $retval;
}

=item $gene->createGroup($groupID[, $filter]);

B<Description:> Add a pre-existing feature to specified group identified by 
    group id.  If the group exists, the feature(s) will be added to it.

B<Parameters:> $groupID - Identifier of the group
               $filter - see format in Gene.pm perldoc (hash ref)

B<Returns:> The number of features added to the specified group.

=cut

sub createGroup {
    #Creates a gene group
    my ($self, $groupID, $filter) = @_;  
    $filter = {'all' => 1} unless($filter);
    return $self->addGenesToGroup($groupID, $filter);
}

=item $gene->removeFromGroup($groupID[, $filter]);

B<Description:> Removes features from specified group based on the filter.
    If no filter is provided, the group is removed.

B<Parameters:> $groupID - Identifier of the group
               $filter - see Gene.pm perldoc for format

B<Returns:> The number of features deleted.

=cut
sub removeFromGroup {
    #Removes somemthing from group
    my ($self, $groupID, $filter) = @_;
    die("$groupID is not a valid group id\n") unless(exists($self->{'groups'}->{$groupID}));
    my $delete = $self->filterFunction($filter) if($filter);
    my $retval = scalar @{$delete};
    
    unless($filter) {
        $retval = scalar (keys %{$self->{'groups'}->{$groupID}});
        delete($self->{'groups'}->{$groupID});
    } else {
        foreach my $feat (@{$delete}) {
            my $id2Del = $feat->id;
            delete($self->{'groups'}->{$groupID}->{$id2Del}) or
                $retval--;
        }
    }

    return $retval;

}

=item $gene->getGroup($groupId);

B<Description:> To retrieve all the features in the group.

B<Parameters:> $groupId - Unique identifier of a group

B<Returns:> Nothing if group does not exist, a hash ref containing all
    the features if it does

=cut
sub getGroup {
    my ($self, $groupId) = @_;
    return [values %{$self->{'groups'}->{$groupId}}];    
}


######################Printing Options################################

=item $gene->printGroup($groupID);

B<Description:> Prints all members of a specified group

B<Parameters:> $groupID - Identifier for the group to be printed.

B<Returns:> The number of features printed

=cut

sub printGroup {
    my ($self, $groupID) = @_;
    &printFeatures([values %{$self->{'groups'}->{$groupID}}]);
}


=item printFeatures($feats);

B<Description:> Class subroutine that will print a set of features

B<Parameters:> $feats - hashRef (for format see the return value of an accessor 
                                 subroutine (getWhatever)).

B<Returns:> The number of features printed.

=cut

sub printFeatures {
    my $feats = shift;
    my $retval = scalar @{$feats};
    my ($k, $v);
    foreach $v (@{$feats}) {
        print "ID     :: ".$v->id."\n";
        print "start  :: ".$v->start."\n";
        print "stop   :: ".$v->stop."\n";
        print "strand :: ".$v->strand."\n";
        print "type   :: ".$v->type."\n\n";
        
    }

    return $retval;
}

=item $gene->filterFunction($filter);

B<Description:> Will filter a genes features based on the filter hash reference.
   For a description of the format of this filter please see the documentation of
   this module.

B<Parameters:> $filter - hashreference

B<Returns:> A hash of features based on the filter

=cut

sub filterFunction {
    my ($self, $filter, $checkThese) = @_;
    $checkThese = $self->{'features'} unless($checkThese);
    my $retval = [];
    
    return $checkThese if($filter->{'all'} || ($filter->{'and'} && scalar(keys %{$filter}) == 1));
    my $and = 1 if($filter->{'and'});

    my ($fTerm, $fVal);

    while(($fTerm,$fVal) = each(%{$filter})) {

        next if($fTerm !~ /(id|type|stop|start|strand)/);
        foreach my $feature (@{$checkThese}) {
            if($fTerm eq 'stop' || $fTerm eq 'start') {
                my ($sign, $val) = ($1,$2) if($filter->{$fTerm} =~ /(.*?)(\d+)/);
                return [] unless($val);

                #Okay, so this craziness right here.  If the sign is negative (meaning the value of the feature
                #should be less than the filter value (ex '-2' means 'less than two')), we should get a postive 
                #result if we subtract the filter value from the feature value.  And opposite if it's not bigger.
                #It should be exactly the opposite if the sign in +.  So by changing the values of a and b, we 
                #can still subtract one from the other and get the correct value.
                my ($a, $b) = ($sign eq '-') ? ($val, $feature->$fTerm()) : ($feature->$fTerm(), $val);
                push(@{$retval}, $feature) if($a-$b > 0 && $sign);
                push(@{$retval}, $feature) if($a-$b == 0 and !($sign));
                    
                
            } else {
                push(@{$retval}, $feature) if($feature->$fTerm() eq $filter->{$fTerm});
            }

            
        }
        
        delete($filter->{$fTerm});
        $retval = $self->filterFunction($filter, $retval) if($and);
        
    }

    return  $retval;
    
}


############################# "PRIVATE" Sub Routines ###################
sub _equal {
    my ($one, $two);
    my $retval = 1;
    
    return $retval;
    
}

sub _init {
    my $self = shift;
    my @params = @_;
    
    if(@params > 0) {
        ($self->{'id'},
         $self->{'start'},
         $self->{'stop'},
         $self->{'strand'},
         $self->{'seq'},
         $self->{'score'}) = @params;       
    }
    
}

sub _mergHashRef {
    my ($one, $two) = @_;

    my ($k, $v);
    
    while(($k, $v) = each(%{$one})) {
        $two->{$k} = $v unless(exists($two->{$k}));
    }

    return $two;
}
1;
