#! /local/perl/bin/perl

use strict;
use warnings;
use lib "/export/CVS/bsml/src"; 
use Pod::Usage;
use BSML::BsmlParserSerialSearch;
use BSML::BsmlReader;
use BSML::BsmlBuilder;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use Data::Dumper;

my %options = ();
my $results = GetOptions( \%options, 'bsmlSearchFile|b=s', 'bsmlMapDir|m=s', 'help|h', 'man' ) || pod2usage();

my $mapping_Parser = new BSML::BsmlParserSerialSearch( SequenceCallBack => \&sequenceHandler );
my $reader = new BSML::BsmlReader;

my $lengthLookUp = {};
my $sequenceLookup = {};
my $offsetLookup = {};

# Process the mapping files, building lookups for sequence length, sequence references, and sequence offsets

print "$options{'bsmlMapDir'}\n";

foreach my $bsmlMap ( <$options{'bsmlMapDir'}/*.bsml> )
{
    $mapping_Parser->parse( $bsmlMap );
}

sub mapping_sequenceHandler
{
    my $seq = shift;
    
    # if the sequence has a numbering it's a chunk, else it's a refernce sequence

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

my $search_Parser = new BSML::BsmlParserSerialSearch( SequenceCallBack => \&search_sequenceHandler, AlignmentCallBack => \&search_alignmentHandler );

sub search_sequenceHandler
{
    my $seqref = shift;

    
}

sub search_alignmentHandler
{
    my $alnref = shift;
}
