#! /local/perl/bin/perl
#

use strict;
use warnings;
use Pod::Usage;
use BSML::BsmlParserSerialSearch;
use BSML::BsmlReader;
use BSML::BsmlBuilder;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use Data::Dumper;

my %options = ();
my $results = GetOptions( \%options, 'bsmlMapSearchFile|b=s', 'bsmlSearchDir|d=s', 'bsmlMapDir|m=s', 'bsmlOutputDir|o=s', 'help|h', 'man' ) || pod2usage();

my $mapping_Parser = new BSML::BsmlParserSerialSearch( SequenceCallBack => \&mapping_sequenceHandler );
my $post_mapping_Parser = new BSML::BsmlParserSerialSearch( SequenceCallBack => \&post_mapping_sequenceHandler );

my $reader = new BSML::BsmlReader;

my $lengthLookUp = {};
my $sequenceLookup = {};
my $offsetLookup = {};

my $currentAssembly = '';
my $currentChunkList = [];

my $bsmlDoc;
my $analysis_written = 0;

# Process the mapping files, building lookups for sequence length, sequence references, and sequence offsets

foreach my $bsmlMap ( <$options{'bsmlMapDir'}/*.bsml> )
{
    $mapping_Parser->parse( $bsmlMap );
}

sub mapping_sequenceHandler
{
    my $seq = shift;
    
    # if the sequence has a numbering it's a chunk, else it's a reference sequence

    if( my $numbering = $reader->readNumbering( $seq ) )
    {
	my $id = $seq->returnattr( 'id' );
	my $length = $seq->returnattr( 'length' );
	
	my $seqref = $numbering->{'seqref'};
	my $offset = $numbering->{'refnum'};

	$sequenceLookup->{$id} = $seqref;
	$offsetLookup->{$id} = $offset;
    }
    else
    {
	my $id =  $seq->returnattr( 'id' );
	my $length = $seq->returnattr( 'length' );
	
	$lengthLookUp->{$id} = $length;
    }
}

sub post_mapping_sequenceHandler
{
  my $seq = shift;
    
    # if the sequence has a numbering it's a chunk, else it's a reference sequence

    if( my $numbering = $reader->readNumbering( $seq ) )
    {
	push( @{$currentChunkList}, $seq->returnattr( 'id' ));
    }
    else
    {
	$currentAssembly =  $seq->returnattr( 'id' );
    }
}

my $search_Parser = new BSML::BsmlParserSerialSearch( SequenceCallBack => \&search_sequenceHandler, AlignmentCallBack => \&search_alignmentHandler, AnalysisCallBack => \&analysis_Handler );

if( !( $options{'bsmlMapSearchFile'} ) )
{
    foreach my $bsmlMap ( <$options{'bsmlMapDir'}/*.bsml> )
    {
	$currentAssembly = '';
	$currentChunkList = [];
	
	$post_mapping_Parser->parse( $bsmlMap );

	$bsmlDoc = new BSML::BsmlBuilder;
	$bsmlDoc->makeCurrentDocument();
	
	my $search = '';
	
	foreach my $chunk ( <$options{'bsmlSearchDir'}/${currentAssembly}_*>)
	{
	    $chunk =~ /$options{'bsmlSearchDir'}\/${currentAssembly}_([\d]*)_([\S]*).bsml/;
	    $search = $2;

	    $search_Parser->parse( $chunk );
	}

	$bsmlDoc->write( "$options{'bsmlOutputDir'}/${currentAssembly}_${search}.bsml" );
    }
}
else
{
   
    my $bsmlMap = $options{'bsmlMapSearchFile'};

    $currentAssembly = '';
    $currentChunkList = [];
	
    $post_mapping_Parser->parse( $bsmlMap );
    
    $bsmlDoc = new BSML::BsmlBuilder;
    $bsmlDoc->makeCurrentDocument();
    
    my $search = '';
	
    foreach my $chunk ( <$options{'bsmlSearchDir'}/${currentAssembly}_*>)
    {
	print "$chunk\n";
	$chunk =~ /$options{'bsmlSearchDir'}\/${currentAssembly}_([\d]*)_([\S]*).bsml/;
	$search = $2;
	
	$search_Parser->parse( $chunk );
    }
    
    $bsmlDoc->write( "$options{'bsmlOutputDir'}/${currentAssembly}_${search}.bsml" );
}

sub search_sequenceHandler
{
    my $seqref = shift;

    # For each sequence, look up the identifier in the map. If the id maps to a super-sequence, check/add the supersequence to the reAssembled document. If the 
    # id is not in the map, pass it through.


    my $seq_length = $seqref->returnattr( 'length' );
    my $seq_title = $seqref->returnattr( 'title' );
    my $seq_molecule = $seqref->returnattr( 'molecule' );
    my $seq_id = $seqref->returnattr( 'id' );

    if( my $mapped_id = $sequenceLookup->{$seq_id} )
    {
	$seq_id = $mapped_id;
	$seq_length = $lengthLookUp->{$mapped_id};
    }

    my $assmblId = $seqref->returnBsmlAttr( 'ASSEMBLY' );

    if( my $mapped_id = $sequenceLookup->{$assmblId} )
    {
	$assmblId = $mapped_id;
    }

    if( !( BSML::BsmlDoc::BsmlReturnDocumentLookup( $seq_id ) ) )
    {
	my $newSeq = $bsmlDoc->createAndAddExtendedSequenceN( id => $seq_id,
							  title => $seq_id, 
							  molecule => $seq_molecule,
							  length => $seq_length );

	$newSeq->addBsmlAttr( 'ASSEMBLY', $assmblId );
    }
}

sub search_alignmentHandler
{
    my $alnref = shift;
    my $rhash = $reader->readSeqPairAlignment( $alnref );

    my $compseq = $rhash->{'compseq'};
    my $refseq = $rhash->{'refseq'};
    my $refoffset = 0;
    my $compoffset = 0;

    if( my $mapped_refseq = $sequenceLookup->{$refseq} )
    {
	$refseq = $mapped_refseq;
	$refoffset = $offsetLookup->{$refseq};
	$rhash->{'reflength'} = $lengthLookUp->{$refseq};
    }

    if( my $mapped_compseq = $sequenceLookup->{$compseq} )
    {
	$compseq = $mapped_compseq;
	$compoffset = $offsetLookup->{$compseq};
	$rhash->{'complength'} = $lengthLookUp->{$compseq};
    }
    
     my $aln = $bsmlDoc->createAndAddSequencePairAlignment( 'refseq' => $refseq,
							'compseq' => $compseq,
							'reflength' => $rhash->{'reflength'},
							'complength' => $rhash->{'complength'}
							);

    foreach my $alnrun ( @{$rhash->{'seqPairRuns'}} )
    {
	my $run = $bsmlDoc->createAndAddSequencePairRun( 'alignment_pair' => $aln,
							 'refpos' => $refoffset + $alnrun->{'refpos'},
							 'comppos' => $compoffset + $alnrun->{'comppos'},
							 'runlength' => $alnrun->{'runlength'},
							 'refcomplement' => $alnrun->{'refcomplement'},
							 'comprunlength' => $alnrun->{'comprunlength'},
							 'compcomplement' => $alnrun->{'compcomplement'},
							 'runscore' => $alnrun->{'runscore'},
							 'runprob' => $alnrun->{'runprob'} );

	foreach my $qual (keys(%{$alnrun}))
	{
	    if( $qual ne 'refpos' &&
		$qual ne 'comppos' &&
		$qual ne 'runlength' &&
		$qual ne 'refcomplement' &&
		$qual ne 'comprunlength' &&
		$qual ne 'compcomplement' &&
		$qual ne 'runscore' &&
		$qual ne 'runprob' &&
		$qual ne 'bsmlRef' )
	    {
		$run->addBsmlAttr( $qual, $alnrun->{$qual} );
	    }
	}
    }
}



sub analysis_Handler
{
    my $analysis = shift;
    my $rhash = $reader->readElement( $analysis );

    
    if( $analysis_written == 0 )
    {
	my $new_analysis = $bsmlDoc->createAndAddAnalysis();

	foreach my $attr (keys(%{$rhash}))
	{
	    $new_analysis->addBsmlAttr( $attr, $rhash->{$attr} );
	}

	$analysis_written = 1;
    }
}
