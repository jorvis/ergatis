#! /usr/local/bin/perl

use strict;
use warnings;

use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);

use lib '/export/CVS/bsml/src';
use BSML::BsmlBuilder;

my %options = ();
my $results = GetOptions( \%options, 'tilingPath|t=s', 'outFile|o=s' );

open( TILINGS, $options{'tilingPath'} ) or die "Unable to open $options{'tilingPath'}\n";

my $bsmlDoc = new BSML::BsmlBuilder;

# add the reference sequence to the Tiling Document

my $line = <TILINGS>;

if( !($line) )
{
    die "Tiling path contains no reference sequence\n";
}
else
{
    # Grab the refseq id from the first line
    $line =~ />([\S]*)[\s]([\d]*)/;

    # add the reference sequence to the bsml tiling
    my $newSeq = $bsmlDoc->createAndAddExtendedSequenceN( id => $1,
							  title => $1,
							  molecule => 'dna',
							  length => $2 );

    $newSeq->addBsmlAttr( 'TYPE', 'REFERENCE_ASSEMBLY' );
}

# add the tilingpath seq to the Tiling Document

my $newSeq = $bsmlDoc->createAndAddExtendedSequenceN( id => 'tiling',
						      title => 'tiling',
						      molecule => 'dna',
						      length => '-1' );


my $rank = 0;

# add each tiled sequence to the Tiling Document

while( my $line = <TILINGS> )
{
    my @tile = split( "\t", $line );
    chomp( $tile[7] );

    my $assmbl_id = $tile[7];
    my $orientation = $tile[6];
    my $ascending = 0;

    if( $orientation eq '+' )
    {
	$ascending = 1;
    }
    else
    {
	$ascending = 0;
    }

    my $newSeq = $bsmlDoc->createAndAddExtendedSequenceN( id => $assmbl_id,
							  title => $assmbl_id,
							  molecule => 'dna',
							  length => $tile[3] );

    $newSeq->addBsmlAttr( 'TYPE', 'TILED_ASSEMBLY' );

    $bsmlDoc->createAndAddNumbering( seq => $newSeq,
				     seqref => 'tiling',
				     refnum => $tile[0],
				     ascending => $ascending );
    
}

$bsmlDoc->write( $options{'outFile'} );
