#! /usr/local/bin/perl

# use lib "/export/CVS/bsml/src";
use strict;
use warnings;
BEGIN {
use BSML::BsmlBuilder;
}
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);

my %options = ();
my $results = GetOptions( \%options, 'mumsAlignFile|m=s', 'outFile|o=s', 'ref_assmbl|r=s', 'help|h', 'man' );

my $mumsAlignFile = $options{'mumsAlignFile'};
my $refAssmbl = $options{'ref_assmbl'};
my $outputName = $options{'outFile'};

my $bsmlDoc = new BSML::BsmlBuilder();

open( ALIGN, $mumsAlignFile ) or die "Unable to open $mumsAlignFile\n" ;

my $assembly_mums = []; 
my $assemblyId = '';
my $revcom = 0;
my $revLength;

while( my $line = <ALIGN> )
{

    if( $line =~ /^> ([\S]*)/ )
    {
	if( $assemblyId && $assembly_mums->[0])
	{

	    my $arr_length = @{$assembly_mums};

	    for( my $i=0; $i<$arr_length; $i++ ) 
	    {
		my $mumline = $assembly_mums->[$i];

		### (INSERTION/DELETION)
       
		if ($mumline->[4] =~ /-/ )
		{
		    my $ref_del_start = $mumline->[2];
		    my $query_del_start = $mumline->[3];
		    my $insertion = $mumline->[5];


		    LINE: while( ++$i < $arr_length )
		    {
			$mumline = $assembly_mums->[$i];

			if( !($mumline->[2] == $ref_del_start))
			{
			    $mumline = $assembly_mums->[--$i];

			    last LINE;
			}
		    
			$insertion .= $mumline->[5];
		    }

		    my $ref_del_end = $mumline->[2];
		    my $query_del_end = $mumline->[3];


		    # Space Based coordinate transformation

		    my $refpos_space_coord;
		    my $querypos_start_space_coord;
		    my $query_length_space_coord;

		    if( $revcom == 1 )
		    {
			$refpos_space_coord = $ref_del_start - 1;
			$querypos_start_space_coord = $revLength - $query_del_end;
			$query_length_space_coord = ( $query_del_end - $query_del_start ) + 1;
		    }
		    else
		    {
			$refpos_space_coord = $ref_del_start - 1;
			$querypos_start_space_coord = $query_del_start - 1;
			$query_length_space_coord = ( $query_del_end - $query_del_start ) + 1;
		    }


		    my $aln = $bsmlDoc->createAndAddSequencePairAlignment(
									  refseq => $refAssmbl,
									  compseq => $assemblyId,
									  class => 'refdel',
									  force => 1 );

		    my $run = $bsmlDoc->createAndAddSequencePairRun(
								    alignment_pair => $aln,
								    refpos => $refpos_space_coord,
								    runlength => 0,
								    refcomplement => 0,
								    comppos => $querypos_start_space_coord,
								    comprunlength => $query_length_space_coord,
								    compcomplement => $revcom );


		    my $refseq = BSML::BsmlDoc::BsmlReturnDocumentLookup( $refAssmbl );
		    $refseq->addattr( 'class', 'assembly' );
		    $refseq->addattr( 'molecule', 'dna' );

		    my $compseq = BSML::BsmlDoc::BsmlReturnDocumentLookup( $assemblyId );
		    $compseq->addattr( 'class', 'assembly' );
		    $compseq->addattr( 'molecule', 'dna' );
		    

		    $aln->addattr( 'class', 'deletion' );
		    $run->addattr( 'alignment', "-..$insertion" );
		}
       
		if ($mumline->[5] =~ /-/ )
		{
		    my $ref_ins_start = $mumline->[2];
		    my $query_ins_start = $mumline->[3];

		    my $insertion = $mumline->[4];

		    LINE: while( ++$i < $arr_length )
		    {
			$mumline = $assembly_mums->[$i];

			if( !($mumline->[3] == $query_ins_start))
			{
			    $mumline = $assembly_mums->[--$i];

			    last LINE;
			}
		    
			$insertion .= $mumline->[4];
		    }

		    my $ref_ins_end = $mumline->[2];
		    my $query_ins_end = $mumline->[3];



		    my $refpos_space_coord;
		    my $querypos_start_space_coord;
		    my $ref_length_space_coord;

		    if( $revcom == 1 )
		    {
			$refpos_space_coord = $ref_ins_start - 1;
			$ref_length_space_coord = ( $ref_ins_end - $ref_ins_start ) + 1;
			$querypos_start_space_coord = $revLength - $query_ins_end + 1;
		    }
		    else
		    {
			$refpos_space_coord = $ref_ins_start - 1;
			$querypos_start_space_coord = $query_ins_end - 1;
			$ref_length_space_coord = ( $ref_ins_end - $ref_ins_start ) + 1;
		    }


		    my $aln = $bsmlDoc->createAndAddSequencePairAlignment(
									  refseq => $refAssmbl,
									  compseq => $assemblyId,
									  class => 'refins',
									  force => 1 );

		    my $run = $bsmlDoc->createAndAddSequencePairRun(
								    alignment_pair => $aln,
								    refpos => $refpos_space_coord,
								    refcomplement => 0,
								    runlength => $ref_length_space_coord,
								    comppos => $querypos_start_space_coord,
								    comprunlength => 0,
								    compcomplement => $revcom );

		    my $refseq = BSML::BsmlDoc::BsmlReturnDocumentLookup( $refAssmbl );
		    $refseq->addattr( 'class', 'assembly' );
		    $refseq->addattr( 'molecule', 'dna' );
		    
		    my $compseq = BSML::BsmlDoc::BsmlReturnDocumentLookup( $assemblyId );
		    $compseq->addattr( 'class', 'assembly' );
		    $compseq->addattr( 'molecule', 'dna' );

		    $aln->addattr( 'class', 'insertion' );
		    $run->addattr( 'alignment', "$insertion..-" );
		}
		
		if ($mumline->[4] =~ /[agct]/ && $mumline->[5] =~ /[agct]/ )
		{	
		    my $refPos = $mumline->[2];
		    my $asblPos = $mumline->[3];

		    my $refBase = $mumline->[4];
		    my $asblBase = $mumline->[5];


		    my $refpos_space_coord;
		    my $querypos_space_coord;

		    if( $revcom == 1 )
		    {
			$refpos_space_coord = $refPos - 1;
			$querypos_space_coord = $revLength - $asblPos;
		    }
		    else
		    {
			$refpos_space_coord = $refPos - 1;
			$querypos_space_coord = $asblPos - 1;
		    }

		    my $aln = $bsmlDoc->createAndAddSequencePairAlignment(
									  refseq => $refAssmbl,
									  compseq => $assemblyId,
									  class => 'snp',
									  force => 1 );

		    my $run = $bsmlDoc->createAndAddSequencePairRun(
								    alignment_pair => $aln,
								    refpos => $refpos_space_coord,
								    runlength => 1,
								    refcomplement => 0,
								    comppos => $querypos_space_coord,
								    comprunlength => 1,
								    compcomplement => $revcom );


		    my $refseq = BSML::BsmlDoc::BsmlReturnDocumentLookup( $refAssmbl );
		    $refseq->addattr( 'class', 'assembly' );
		    $refseq->addattr( 'molecule', 'dna' );
		    
		    my $compseq = BSML::BsmlDoc::BsmlReturnDocumentLookup( $assemblyId );
		    $compseq->addattr( 'class', 'assembly' );
		    $compseq->addattr( 'molecule', 'dna' );

		    $aln->addattr( 'class', 'SNP' );
		    $run->addattr( 'alignment', "$refBase..$asblBase" );
		} 
	    }

	    $assembly_mums = [];
	}

	$assemblyId = $1;

	if( $assemblyId =~ /revcom:([\d]*)/ )
	{
	    $revcom = 1;
	    $revLength = $1;
	    $assemblyId =~ s/revcom:([\d]*)//;
	}
	else
	{
	    $revcom = 0;
	}

	$assemblyId =~ s/ID//;

	if( $line =~ /Reverse/ )
	{
	    $assemblyId = 'SKIP';
	}
	
	next;
    }

    if( !($assemblyId eq 'SKIP') )
    {
	chomp( $line );
	my @mum = split( " ", $line );
	my $mumref = \@mum;

	# Verify that the six columns have been populated 
	
	if( ( $mum[0] && $mum[1] && $mum[2] && $mum[3] && $mum[4] && $mum[5] ))
	{
	   push( @{$assembly_mums}, $mumref );
	}
	      
	
    }

	

}

$bsmlDoc->write( $outputName );
