#!/usr/bin/perl
BEGIN{foreach (@INC) {s/\/usr\/local\/packages/\/local\/platform/}};
use lib (@INC,$ENV{"PERL_MOD_DIR"});
no lib "$ENV{PERL_MOD_DIR}/i686-linux";
no lib ".";

use strict;
use BSML::BsmlParserSerialSearch;

use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);

my %options = ();
my $results = GetOptions( \%options, 'cogFile|c=s', 'bsmlModelList|m=s', 'outputDir|o=s', 'maxCogSeqCount|s=s', 'extension|e=s');

my $cogFile = $options{'cogFile'};
my $bsmlModelList = $options{'bsmlModelList'};
my $outDir = $options{'outputDir'};
my $maxCogSeqCount = $options{'maxCogSeqCount'};
if($options{'extension'} eq ''){
    $options{'extension'} = 'fsa';
}

my $Prot = {};

# set up a serial parser to parse sequence elements, bypassing feature tables for 
# efficiency.

my $seqParser = new BSML::BsmlParserSerialSearch( ReadFeatureTables => 0, SequenceCallBack => \&createPolypeptideLookup);

#Get rid of trailing slashes in directory names

$bsmlModelList =~ s/\/+$//;
$outDir =~ s/\/+$//; 

if( !$bsmlModelList )
{
    die "no Bsml Directory specified\n";
}
else
{
    if( ! -e $bsmlModelList )
    {
	die "could not open list file: $bsmlModelList\n";
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

open FILE, $bsmlModelList or die "Can't open file $bsmlModelList";
while(my $bsmlFile=<FILE>)
{
    chomp $bsmlFile;
    if(-e $bsmlFile){
	print STDOUT "Parsing $bsmlFile\n";
	$seqParser->parse( $bsmlFile );
    }
}

open( INPUTCOGS, "<$cogFile" ) or die "could not open $cogFile.\n";

my $cog = '';
my $list = [];

while( my $line = <INPUTCOGS> )
{
    if( $line =~ /^\t([\S]*)/ )
    {
	# A new sequence has been found and added to the current cog.
	print STDERR "Found line $line\n";
	push( @{$list}, $1 );
    }

    if( $line =~ /^COG\s+=\s+([^,\s]+)/ )
    {
	# A new cog has been encountered, flush the previous to disk if present
	print STDERR "Found cog $line $cog\n";
	if( $cog){
	    print STDERR "Outputing cog $cog\n";
	    if(scalar(@{$list})>1){
		if(@{$list} <= $maxCogSeqCount){
		    open( OUTFILE, ">$outDir/$cog.$$.$options{'extension'}" ) or die "could not open $cog.$options{'extension'}\n";
		    foreach my $seq ( @{$list} )
		    {
			print OUTFILE ">$seq\n";
			print OUTFILE $Prot->{$seq}."\n";
		    }
		    close( OUTFILE );
		}
		else{
		    open( OUTFILE, ">$outDir/$cog.$$.$options{'extension'}" ) or die "could not open $cog.$options{'extension'}\n";
		    foreach my $seq ( @{$list} )
		    {
			print OUTFILE ">$seq\n";
			print OUTFILE "X\n";
		    }
		    close( OUTFILE );
		}
	    }
	}

	$cog = $1;
	$list = [];
    }
}
if( $cog){
    print STDERR "Outputing cog $cog\n";
    if(scalar(@{$list})>1){
	if(@{$list} <= $maxCogSeqCount){
	    open( OUTFILE, ">$outDir/$cog.$$.$options{'extension'}" ) or die "could not open $cog.$options{'extension'}\n";
	    foreach my $seq ( @{$list} )
	    {
		print OUTFILE ">$seq\n";
			print OUTFILE $Prot->{$seq}."\n";
	    }
	    close( OUTFILE );
	}
	else{
	    open( OUTFILE, ">$outDir/$cog.$$.$options{'extension'}" ) or die "could not open $cog.$options{'extension'}\n";
	    foreach my $seq ( @{$list} )
	    {
		print OUTFILE ">$seq\n";
		print OUTFILE "A\n";
	    }
	    close( OUTFILE );
	}
    }
}

sub createPolypeptideLookup
{
    my $seqRef = shift;

    # We're only interested in polypeptide sequences for all-vs-all and pblast

    if( ($seqRef->returnattr( 'molecule' ) eq 'aa') || ($seqRef->returnattr('class') eq 'polypeptide'))
    {
	my $seq = $seqRef->subSequence(-1,0,0);
	
	my $identifier;
	if(defined $seqRef->{'BsmlSeqDataImport'}){
	    $identifier = $seqRef->{'BsmlSeqDataImport'}->{'identifier'};
	}
	else{
	    $identifier = $seqRef->returnattr( 'id' );
	}

	if( $identifier && $seq )
	{
	    $Prot->{$identifier} = $seq;
	}
    }
}







