#! /local/perl/bin/perl

use lib '/export/CVS/bsml/src';
use strict;
use warnings;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use BSML::BsmlBuilder;
use Data::Dumper;

################
#
# Encode MSF (Clustal Multiple Sequence Format) files into BSML
#
#
# Arguments 
#   msfdir - directory of clustal files
#   outputFile - bsml output file
#
#
#

my $builder = new BSML::BsmlBuilder;

my %options = ();
my $results = GetOptions( \%options, 'msfdir|d=s', 'outputFile|o=s' );

my $msfDir = $options{'msfdir'};
my $outFile = $options{'outputFile'};

$msfDir =~ s/\/+$//; #remove trailing slash if present

foreach my $msfFile (<$msfDir/*.clustal>)
{
    print "  - $msfFile\n";
    open( INFILE, "<$msfFile" ) or die "Could not open $msfFile\n";

    my $msfLength;
    my $msfType;
    my $seqs = {};
    my $count = 1;

    while( my $line = <INFILE> )
    {
	if( $line =~ /MSF:[\s]*([\d]*)[\s]*Type:[\s]*([PN])[\s]*Check:[\s]*[\d]*[\s]*../ )
	{
	    $msfLength = $1;
	    
	    $msfType = 'protein' if( $2 eq 'P' );
	    $msfType = 'nucleotide' if( $2 eq 'N' );
	}

	if( $line =~ /Name:[\s]*([\S]*)[\s]*oo[\s]*Len:[\s]*([\d]*)[\s]*Check:[\s]*[\d]*[\s]*Weight:[\s]*([\d.]*)\n/ ) 
	{
	    $seqs->{$1}->{'length'} = $2;
	    $seqs->{$1}->{'weight'} = $3;
	    $seqs->{$1}->{'alignment'} = [];
	    $seqs->{$1}->{'num'} = $count++;
	} 

	foreach my $name ( keys( %{$seqs} ) )
	{
	    if( $line =~ /^$name/ )
	    {
		push( @{$seqs->{$name}->{'alignment'}}, $line );
	    }
	}
    }

    # Prevent encoding of empty objects if empty clustal files are encountered
    if( keys(%{$seqs}) >= 1 )
	{
	    my $table = $builder->createAndAddMultipleAlignmentTable( 'molecule-type' => $msfType );
	    my $summary = $builder->createAndAddAlignmentSummary( 'multipleAlignmentTable' => $table,
								  'seq-type' => $msfType,
								  'seq-format' => 'msf' );

	    my $aln = $builder->createAndAddSequenceAlignment( 'multipleAlignmentTable' => $table );
	    my $sequences_tag = '';
	    
	    foreach my $seq (keys( %{$seqs} ))
	    {
		my $seqnum = $seqs->{$seq}->{'num'};
		my $length = $seqs->{$seq}->{'length'};
		my $name = $seq;
    
		$builder->createAndAddAlignedSequence( alignmentSummary => $summary,
						       seqnum => $seqnum,
						       length => $length,
						       name => $name );

   
		
		$builder->createAndAddSequenceData( 'sequenceAlignment' => $aln,
						    'seq-name' => $name,
						    'seq-data' => join( '', @{ $seqs->{$seq}->{'alignment'} } ));

		$sequences_tag .= "$seqnum:";
	    }

	    chop( $sequences_tag );
	    
	    $aln->addattr( 'sequences', $sequences_tag );
	}
}

	$builder->write( $outFile );

