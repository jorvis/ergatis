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

B<--outFile, -o>        [REQUIRED] output BTAB file

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

use strict;
use warnings;
use Pod::Usage;
use BSML::BsmlParserSerialSearch;
use BSML::BsmlReader;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use Data::Dumper;

my $pairsCount = 0;

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
 
my %options = ();
my $results = GetOptions( \%options, 'bsmlSearchDir|b=s', 'bsmlModelDir|m=s', 'bsmlJaccardDir|j=s', 'jaccardFilter|f=s', 'outFile|o=s', 'pvalcut|p=s', 'help|h', 'man' ) || pod2usage();

my $bsmlSearchDir = $options{'bsmlSearchDir'};
my $bsmlModelDir = $options{'bsmlModelDir'};my $outFile = $options{'outFile'};
my $PVALCUTOFF = $options{'pvalcut'};

# associative array to translate cds identifiers to protein ids.
my $cds2Prot = {};

# The alignment parser handles the pairwise alignments encoded in the search directory. The feature parser
# creates a lookup table mapping protein sequence identifiers to their genome. 

my $alnParser = new BSML::BsmlParserSerialSearch( AlignmentCallBack => \&alignmentHandler );
my $featParser = new BSML::BsmlParserSerialSearch( FeatureCallBack => \&featureHandler, GenomeCallBack => \&genomeHandler );

if( exists($options{'help'}) || exists($options{'man'}) )
{
    pod2usage({-exitval => 1, -verbose => 2, -output => \*STDOUT});
}

if( !$PVALCUTOFF )
{
    die "ERROR: No pvalue cutoff has been specified";
}

# If Jaccard data has been specified, use it for equivalence class filtering

my $jaccardClusterHash = {};  #Associate Jaccard clusters by sequence id
my $jaccardRepSeqHash = {};   #Associate a representative sequence for each cluster
my $jaccardClusterCount = 0;

if( $options{'bsmlJaccardDir'} && $options{'jaccardFilter'} == 1 )
{
    print "Filtering\n";

    my $multiAlnParser = new BSML::BsmlParserSerialSearch( MultipleAlignmentCallBack => \&multipleAlignmentHandler );
    foreach my $bsmlFile (<$options{'bsmlJaccardDir'}/*.bsml>)
    {
	$multiAlnParser->parse( $bsmlFile );
    }
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

	    $jaccardClusterHash->{$name} = $jaccardClusterCount;
	    if( $seqCount == 0 )
	    {
		$jaccardRepSeqHash->{$jaccardClusterCount} = $name;
	    }
	    $seqCount++;
	}
    }

    $jaccardClusterCount++;
}

if(!$bsmlSearchDir || !$bsmlModelDir )
{
    die "ERROR: BSML directories for search encodings and gene-sequence models not specified\n";
}
else
{
    $bsmlSearchDir =~ s/\/+$//; #remove trailing slash if present
    $bsmlModelDir =~ s/\/+$//; 

    if( ! -d $bsmlSearchDir || ! -d $bsmlModelDir )
    {
	die "could not open directory: $bsmlSearchDir\n";
    }
}



# protein sequence identifer to genome mapping 
my $geneGenomeMap = {};
my $genome = '';

# loop through the documents in the model directory to create the protein genome map

foreach my $bsmlFile (<$bsmlModelDir/*.bsml>)
{
    $featParser->parse( $bsmlFile );
    $genome = '';
}


if(! $outFile )
{
    die "output file not specified\n";
}

open( OUTFILE, ">$outFile" );

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

foreach my $bsmlFile (<$bsmlSearchDir/*.bsml>)
{
    
    # builds the COGS input data structure

    $alnParser->parse( $bsmlFile );

    # print the results

    foreach my $k1 ( keys( %{$COGInput} ) )
    {
	foreach my $k2 (keys( %{$COGInput->{$k1}}))
	{
	    print OUTFILE join("\t", @{$COGInput->{$k1}->{$k2}});
	    print OUTFILE "\n";
	}
    }

    $COGInput = {};
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


    # if compseq is part of a Jaccard equivalence class, only use it if it is the
    # Jaccard reference sequence for the class

    if( defined( my $jId = $jaccardClusterHash->{$compseq} ) )
    {
	if( !($compseq eq $jaccardRepSeqHash->{$jId}) )
	{
	    return;
	}
    }

    # find the highest scoring run for each pairwise alignment

    foreach my $seqPairRun ( @{$aln->returnBsmlSeqPairRunListR} )
    {
	my $runscore = $seqPairRun->returnattr( 'runscore' );
	my $pvalue = $seqPairRun->returnBsmlAttr( 'p_value' );

	if( ($runscore > $bestRunScore) && ($pvalue < $PVALCUTOFF) )
	{
	    $bestRunScore = $runscore;
	    $bestSeqPairRun = $seqPairRun;
	}
    }

    # 
    return if( !($bestSeqPairRun) );

    my $runscore = $bestSeqPairRun->returnattr( 'runscore' );
    my $compGenome = $geneGenomeMap->{$compseq};
    my $refGenome = $geneGenomeMap->{$refseq};

    if( !( $compGenome )) 
    {
	die "$compseq: compseq not found in gene genome map\n";
    }

    if( !( $refGenome ) )
    {
	die "$refGenome: refseq not found in gene genome map\n";
    }

    # alignments to the same genome are not included in COG clusters

    return if( $compGenome eq $refGenome );

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
	    $COGInput->{$refseq}->{$compGenome} = $lref;
	}
    }
    else
    {
	$COGInput->{$refseq}->{$compGenome} = $lref;
    }

    if( $pairsCount > 200 ){ exit(); }
    $pairsCount++;
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

		$cds2Prot->{$cdsId} = $protId;
		$geneGenomeMap->{$protId} = $genome;

		return;
	    }
	}
    }
} 



