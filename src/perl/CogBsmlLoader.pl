#! /local/perl/bin/perl

use lib '/export/CVS/bsml/src';
use strict;
use warnings;
use BSML::BsmlParserSerialSearch;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);

# Preprosess data stored in BSML pairwise alignment documents into BTAB
# structure for COG analysis.

########
# Arguments:
#
# bsmlSearchDir - directory containing BSML pairwise sequence encodings
# bsmlModelDir - directory containing the BSML gene model documents for the search data
# outfile - btab output file
#
#
 
my %options = ();
my $results = GetOptions( \%options, 'bsmlSearchDir|b=s', 'bsmlModelDir|m=s', 'outFile|o=s' );

my $bsmlSearchDir = $options{'bsmlSearchDir'};
my $bsmlModelDir = $options{'bsmlModelDir'};
my $outFile = $options{'outFile'};

my $cds2Prot = {};

# Usage of two parsers is less efficient, but insures that the sequence objects are 
# parsed before the alignments.

my $seqParser = new BSML::BsmlParserSerialSearch( SequenceCallBack => \&createGeneGenomeMap );
my $alnParser = new BSML::BsmlParserSerialSearch( AlignmentCallBack => \&alignmentHandler );
my $featParser = new BSML::BsmlParserSerialSearch( FeatureCallBack => \&featureHandler );

print "Looking for BSML documents in $bsmlSearchDir\n";

$bsmlSearchDir =~ s/\/+$//; #remove trailing slash if present
$bsmlModelDir =~ s/\/+$//; 

if(!$bsmlSearchDir || !$bsmlModelDir )
{
    die "no Bsml Directory specified\n";
}
else
{
    if( ! -d $bsmlSearchDir || ! -d $bsmlModelDir )
    {
	die "could not open directory: $bsmlSearchDir\n";
    }
}

foreach my $bsmlFile (<$bsmlModelDir/*.bsml>)
{
    print "  - $bsmlFile\n";
    $featParser->parse( $bsmlFile );
}


if(! $outFile )
{
    die "output file not specified\n";
}

open( OUTFILE, ">$outFile" );
    
my $geneGenomeMap = {};
my $COGInput = {};

foreach my $bsmlFile (<$bsmlSearchDir/*.bsml>)
{
    $geneGenomeMap = {};
    print "  - $bsmlFile\n";
    $seqParser->parse( $bsmlFile );
    $alnParser->parse( $bsmlFile );

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

sub createGeneGenomeMap
{
    my $seqRef = shift;
    my $gene = $seqRef->returnattr( 'id' );
    my $genome = $seqRef->returnBsmlAttr( 'ASSEMBLY' );
    
    $geneGenomeMap->{$gene} = $genome;
}

sub alignmentHandler
{
    my $aln = shift;
    my $refseq = $aln->returnattr( 'refseq' );
    my $compseq = $aln->returnattr( 'compseq' );
    $refseq = $cds2Prot->{$refseq};
    
    my $bestRunScore = 0;
    my $bestSeqPairRun = undef;

    return if( $compseq eq $refseq );

    foreach my $seqPairRun ( @{$aln->returnBsmlSeqPairRunListR} )
    {
	my $runscore = $seqPairRun->returnattr( 'runscore' );

	if( $runscore > $bestRunScore )
	{
	    $bestRunScore = $runscore;
	    $bestSeqPairRun = $seqPairRun;
	}
    }

    return if( !($bestSeqPairRun) );

    my $runscore = $bestSeqPairRun->returnattr( 'runscore' );
    my $compGenome = $geneGenomeMap->{$compseq};

    if( !( $compGenome )) 
    {
	die "compseq not found in gene genome map\n";
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
	    $COGInput->{$refseq}->{$compGenome} = $lref;
	}
    }
    else
    {
	$COGInput->{$refseq}->{$compGenome} = $lref;
    }
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
		return;
	    }
	}
    }
} 



