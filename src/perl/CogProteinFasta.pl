#! /local/perl/bin/perl

use lib '/export/CVS/bsml/src';
use strict;
use warnings;
use BSML::BsmlParserSerialSearch;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);

my %options = ();
my $results = GetOptions( \%options, 'cogFile|c=s', 'bsmlModelDir|m=s', 'outputDir|o=s', 'maxCogSeqCount|s=s');

my $cogFile = $options{'cogFile'};
my $bsmlModelDir = $options{'bsmlModelDir'};
my $outDir = $options{'outputDir'};
my $maxCogSeqCount = $options{'maxCogSeqCount'};

my $Prot = {};

# set up a serial parser to parse sequence elements, bypassing feature tables for 
# efficiency.

my $seqParser = new BSML::BsmlParserSerialSearch( SequenceCallBack => \&createProteinLookup, ReadBsmlTables => 0 );

#Get rid of trailing slashes in directory names

$bsmlModelDir =~ s/\/+$//;
$outDir =~ s/\/+$//; 

if( !$bsmlModelDir )
{
    die "no Bsml Directory specified\n";
}
else
{
    if( ! -d $bsmlModelDir )
    {
	die "could not open directory: $bsmlModelDir\n";
    }
}

if( !$cogFile )
{
    die "cog file not specified\n";
}

if( !$outDir )
{
    die "output directory not specified\n";
}
else
{
    if( ! -d $outDir )
    {
	mkdir( $outDir );
    }
}


foreach my $bsmlFile (<$bsmlModelDir/*.bsml>)
{
    $seqParser->parse( $bsmlFile );
}

open( INPUTCOGS, "<$cogFile" ) or die "could not open $cogFile.\n";

my $cog = '';
my $list = [];

while( my $line = <INPUTCOGS> )
{
    if( $line =~ /^\t([\S]*)/ )
    {
	# A new sequence has been found and added to the current cog.

	push( @{$list}, $1 );
    }

    if( $line =~ /^COG = ([\d]*),/ )
    {
	# A new cog has been encountered, flush the previous to disk if present

	if( $cog && @{$list} <= $maxCogSeqCount )
	{
	    open( OUTFILE, ">$outDir/$cog.fasta" ) or die "could not open $cog.fasta\n";
	    foreach my $seq ( @{$list} )
	    {
		print OUTFILE ">$seq\n";
		print OUTFILE $Prot->{$seq}."\n";
	    }
	    close( OUTFILE );
	}

	$cog = $1;
	$list = [];
    }
}

sub createProteinLookup
{
    my $seqRef = shift;

    # We're only interested in protein sequences for all-vs-all and pblast

    if( $seqRef->returnattr( 'molecule' ) eq 'aa' )
    {
	my $id = $seqRef->returnattr( 'id' );
	my $seqdat = $seqRef->returnSeqData();

	if( $id && $seqdat )
	{
	    $Prot->{$id} = $seqdat;
	}
    }
}




