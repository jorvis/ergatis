#! /local/perl/bin/perl

=head1  NAME 

CogBsmlLoader.pl  -  Preprosess data stored in BSML pairwise alignment documents into BTAB
structure for COG analysis using best_hits.pl. 

=head1 SYNOPSIS

USAGE:  CogBsmlLoader.pl -m BsmlGeneModelDirectory -b BsmlPairwiseAlignmentDirectory -o OutputFile

=head1 OPTIONS

=over 4

=item *

B<--bsmlModelDir, -m>   [REQUIRED] Dir containing the BSML gene/sequence encodings referenced in the search directory

=item *

B<--bsmlSearchDir, -b>  [REQUIRED] Dir containing the BSML search encodings of pairwise alignments (all_vs_all, blastp)

=item *

B<--output, -o>        [REQUIRED] output BTAB file

=item *

B<--help,-h> This help message

=back

=head1   DESCRIPTION

CogBsmlLoader.pl is designed to preprocess the data contained in a BSML pairwise alignment search 
for COGS analysis. Specifically it identifies the "best hit" per genome for each query gene. 
This data is packaged into the BTAB format for linkage analysis using best_hits.pl  

NOTE:  

Calling the script name with NO flags/options or --help will display the syntax requirement.

=cut

# Preprosess data stored in BSML pairwise alignment documents into BTAB
# structure for COG analysis.

############
# Arguments:
#
# bsmlSearchDir - directory containing BSML pairwise sequence encodings
# bsmlModelDir - directory containing the BSML gene model documents for the search data
# outfile - btab output file
#
#

use strict;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Pod::Usage;
use Workflow::Logger;
use BSML::BsmlParserSerialSearch;
use BSML::BsmlReader;

my %options = ();
my $results = GetOptions( \%options, 
			  'bsmlSearchList|b=s', 
			  'bsmlModelList|m=s', 
			  'bsmlJaccardDir|j=s', 
			  'outfile|o=s', 
			  'pvalcut|p=s', 
			  'log|l=s',
			  'debug=s',
			  'help|h') || pod2usage();

my $logfile = $options{'log'} || Workflow::Logger::get_default_logfilename();
my $logger = new Workflow::Logger('LOG_FILE'=>$logfile,
				  'LOG_LEVEL'=>$options{'debug'});
$logger = Workflow::Logger::get_logger();

# display documentation
if( $options{'help'} ){
    pod2usage( {-exitval=>0, -verbose => 2, -output => \*STDERR} );
}


&check_parameters(\%options);

#MAIN HERE

# associative array to translate cds identifiers to protein ids.
my $cds2Prot = {};

# The alignment parser handles the pairwise alignments encoded in the search directory. The feature parser
# creates a lookup table mapping protein sequence identifiers to their genome. 

my $alnParser = new BSML::BsmlParserSerialSearch( AlignmentCallBack => \&alignmentHandler );
my $featParser = new BSML::BsmlParserSerialSearch( FeatureCallBack => \&featureHandler, GenomeCallBack => \&genomeHandler );

# If Jaccard data has been specified, use it for equivalence class filtering

my $jaccardClusterHash = {};  #Associate Jaccard clusters by sequence id
my $jaccardRepSeqHash = {};   #Associate a representative sequence for each cluster
my $jaccardClusterCount = 0;


$options{'bsmlJaccardDir'} =~ s/\'//g;
$options{'bsmlJaccardDir'} =~ s/\"//g;

if( $options{'bsmlJaccardDir'} && $options{'bsmlJaccardDir'} ne "" )
{
    if(-r $options{'bsmlJaccardDir'}){
	my $multiAlnParser = new BSML::BsmlParserSerialSearch( MultipleAlignmentCallBack => \&multipleAlignmentHandler );
	foreach my $bsmlFile (<$options{'bsmlJaccardDir'}/*.bsml>)
	{
	    $multiAlnParser->parse( $bsmlFile );
	}
    }
    else{
	$logger->logdie("Can't read jaccard dir $options{'bsmlJaccardDir'}");
    }
}


# protein sequence identifer to genome mapping 
my $geneGenomeMap = {};
my $genome = '';

# loop through the documents in the model directory to create the protein genome map

foreach my $bsmlFile (@{&get_list_from_file($options{'bsmlModelList'})}){
    $featParser->parse( $bsmlFile );
    $genome = '';
}


open( OUTFILE, ">$options{'outfile'}" ) or $logger->logdie("Can't open file $options{'outfile'}");

#####################################

# structure for building the COGS input. For each query gene, the COGS analysis expects
# the single best scoring alignment for each reference genome. In BSML terms, COGS expects the
# highest scoring seqpair for each refseq compseq combination where all compseqs 
# are contained in the same genome. 


#  Genome A           Genome B            Genome C
#  compseqA1          compseqB1           compseqC1
#  compseqA2          compseqB2           compseqC2
#  compseqA3                              compseqC3


# If the above represent the sets of reference genes by genome. The following would 
# correspond to the expected output if A1, B2, and C1 were the best hits by genome. 

# refseq -> compseqA1
# refseq -> compseqB2
# refseq -> compseqC1

####################################
    
my $COGInput = {};

foreach my $bsmlFile (@{&get_list_from_file($options{'bsmlSearchList'})}){
    
    # builds the COGS input data structure

    $alnParser->parse( $bsmlFile );

    # print the results

    foreach my $k1 ( keys( %{$COGInput} ) )
    {
	foreach my $k2 (keys( %{$COGInput->{$k1}}))
	{
	    my $member = $COGInput->{$k1}->{$k2}->[0];
	    if(exists $jaccardClusterHash->{$member}){
		$COGInput->{$k1}->{$k2}->[21] = join(',',@{$jaccardRepSeqHash->{$jaccardClusterHash->{$member}}});
	    }
	    print OUTFILE join("\t", @{$COGInput->{$k1}->{$k2}});
	    print OUTFILE "\n";
	}
    }

    $COGInput = {};
}


sub multipleAlignmentHandler
{
    my $aln = shift;
    my $bsml_reader = new BSML::BsmlReader();

    my $maln_ref = $bsml_reader->readMultipleAlignmentTable($aln);

    foreach my $alnSum ( @{ $maln_ref->{'AlignmentSummaries'} } )
    {
	my $seqCount = 0;
	foreach my $alnSeq ( @{ $alnSum->{'AlignedSequences'} } )
	{
	    my $name = $alnSeq->{'name'};
	    $name =~ s/:[\d]*//;

	    # Associate a SeqId with a Jaccard Cluster ID
	    $jaccardClusterHash->{$name} = $jaccardClusterCount;

	    # If this is the first sequence observed in the Jaccard Cluster
	    # identify it as the representative sequence and set it as the 
	    # representative sequence for the Jaccard Cluster ID.
 
	    if( $seqCount == 0 )
	    {
		$jaccardRepSeqHash->{$jaccardClusterCount} = [];
		push @{$jaccardRepSeqHash->{$jaccardClusterCount}},$name;
	    }
	    else{
		push @{$jaccardRepSeqHash->{$jaccardClusterCount}},$name;
	    }
	    $seqCount++;
	}
    }

    $jaccardClusterCount++;
}
sub alignmentHandler
{
    my $aln = shift;
    my $refseq = $aln->returnattr( 'refseq' );
    my $compseq = $aln->returnattr( 'compseq' );

    # if refseq is a CDS identifier, translate it to a protein id. 

    if( my $Trefseq = $cds2Prot->{$refseq} )
    {
	$refseq = $Trefseq;
    }
    
    my $bestRunScore = 0;
    my $bestSeqPairRun = undef;

    # self-self alignments are not included 
    return if( $compseq eq $refseq );

    my $compGenome = $geneGenomeMap->{$compseq};
    my $refGenome = $geneGenomeMap->{$refseq};

    if( !( $compGenome )) 
    {
	$logger->debug("$compseq: compseq not found in gene genome map. skipping.") if($logger->is_debug());
	return;
    }

    if( !( $refGenome ) )
    {
	$logger->debug("$refseq: compseq not found in gene genome map. skipping.") if($logger->is_debug());
	return;
    }

    # alignments to the same genome are not included in COG clusters

    return if( $compGenome eq $refGenome );

     # find the highest scoring run for each pairwise alignment

    foreach my $seqPairRun ( @{$aln->returnBsmlSeqPairRunListR} )
    {
	my $runscore = $seqPairRun->returnattr( 'runscore' );
	my $pvalue = $seqPairRun->returnBsmlAttr( 'p_value' );

	if( ($runscore > $bestRunScore) && ($pvalue < $options{'pvalcut'}) )
	{
	    $logger->debug("Using run with runscore $runscore $pvalue. Previous bestrunscore $bestRunScore. pvalue cutoff $options{'pvalcut'}");
	    $bestRunScore = $runscore;
	    $bestSeqPairRun = $seqPairRun;
	}
    }

    # 
    if( !($bestSeqPairRun) ){
	$logger->warn("Best run not defined");
	return;
    }

    my $runscore = $bestSeqPairRun->returnattr( 'runscore' );
    
    # If compseq (or refseq) is defined in a Jaccard equivalence class identify the class by
    # its reference sequence. 

    if( defined( my $jId = $jaccardClusterHash->{$compseq} ) )
    {
	$logger->debug("Found jaccard cluster $jId for id $compseq. Using $jaccardRepSeqHash->{$jId}->[0] as cluster representative");
	$compseq = $jaccardRepSeqHash->{$jId}->[0];
    }

    if( defined( my $jId = $jaccardClusterHash->{$refseq} ) )
    {
	$logger->debug("Found jaccard cluster $jId for id $refseq. Using $jaccardRepSeqHash->{$jId}->[0] as cluster representative");
	$refseq = $jaccardRepSeqHash->{$jId}->[0];
    }
    
    my $lref = [];
    
    $lref->[0] = $refseq;  #query name
    $lref->[1] = '';       #date
    $lref->[2] = $aln->returnattr( 'reflength' ); #query length
    $lref->[3] = $aln->returnattr( 'method' ); #program
    $lref->[4] = $aln->returnattr( 'compxref' );
    $lref->[5] = $compseq;
    $lref->[6] = $bestSeqPairRun->returnattr( 'refpos' );
    $lref->[7] = $bestSeqPairRun->returnattr( 'refpos' ) + $bestSeqPairRun->returnattr('runlength');
    $lref->[8] = $bestSeqPairRun->returnattr( 'comppos' );
    $lref->[9] = $bestSeqPairRun->returnattr( 'comppos' ) + $bestSeqPairRun->returnattr('comprunlength');
    $lref->[10] = $bestSeqPairRun->returnBsmlAttr( 'percent_identity' );
    $lref->[11] = $bestSeqPairRun->returnBsmlAttr( 'percent_similarity' );
    $lref->[12] = $bestSeqPairRun->returnattr( 'runscore' );
    $lref->[13] = $bestSeqPairRun->returnBsmlAttr( 'chain_number' );
    $lref->[14] = $bestSeqPairRun->returnBsmlAttr( 'segment_number' );
    $lref->[15] = '';
    $lref->[16] = '';
    $lref->[17] = '';
    $lref->[18] = $bestSeqPairRun->returnattr('comprunlength');
    $lref->[19] = $bestSeqPairRun->returnattr('runprob' );
    $lref->[20] = $bestSeqPairRun->returnBsmlAttr( 'p_value' );

   
    if(  $COGInput->{$refseq}->{$compGenome} )
    {
	if(  $COGInput->{$refseq}->{$compGenome}->[12] < $bestRunScore )
	{
	    $logger->debug("$refseq match to $compGenome with score $bestRunScore is highest scoring match.  Previous high score is $COGInput->{$refseq}->{$compGenome}->[12]");
	    $COGInput->{$refseq}->{$compGenome} = $lref;
	}
    }
    else
    {
	$logger->debug("$refseq match to $compGenome is first match found.");
	$COGInput->{$refseq}->{$compGenome} = $lref;
    }
}

sub genomeHandler
{
    my $bsmlGenome = shift;
    my $reader = new BSML::BsmlReader;
    
    my $rhash = $reader->readGenome( $bsmlGenome );

    if( !defined($rhash->{'strain'}) )
    {
	$rhash->{'strain'} = ' ';
    }

    $genome = $rhash->{'genus'}.':'.$rhash->{'species'}.':'.$rhash->{'strain'};
}
	
sub featureHandler
{
    my $feature = shift;

    if( $feature->returnattr('class') eq 'CDS' )
    {
	my $cdsId = $feature->returnattr('id');
	my $protId = '';
	foreach my $link( @{$feature->returnBsmlLinkListR()} )
	{
	    if( $link->{'rel'} eq 'SEQ' )
	    {
		$protId = $link->{'href'};
		$protId =~ s/#//;

		$logger->debug("Adding protein $protId to protein lookup with genome $genome");
		$logger->debug("Adding CDS $cdsId $protId to cds lookup with protiein $protId");

		$cds2Prot->{$cdsId} = $protId;
		$geneGenomeMap->{$protId} = $genome;

		return;
	    }
	}
    }
} 

sub get_list_from_file{
    my($file) = @_;
    my @lines;
    open( FH, $file ) or $logger->logdie("Could not open $file");
    while( my $line = <FH> ){
	chomp($line);
	push @lines, split(',',$line) if($line =~ /\S+/);
    }
    return \@lines;
}

sub check_parameters{
    my ($options) = @_;
    
    if(0){
	pod2usage({-exitval => 2,  -message => "error message", -verbose => 1, -output => \*STDERR});    
    }
}
