#! /local/perl/bin/perl

use strict;
use warnings;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);

my %options = ();
my $results = GetOptions( \%options, 'clustalDir|c=s', 'outputDir|o=s' );

my $clustalDir = $options{'clustalDir'};
my $outputDir = $options{'outputDir'};

$clustalDir =~ s/\/+$//;
$outputDir =~ s/\/+$//; 

die "No Clustal directory specified\n" if( !$clustalDir || ! -d $clustalDir );

if( !$outputDir )
{
    die "No output directory specified\n";
}
else
{
    if( ! -d $outputDir )
    {
	mkdir( $outputDir );
    }
}

foreach my $clustalFile ( <$clustalDir/*.clustal> )
{
    my $bsmlFile = $clustalFile;
    $bsmlFile =~ s/$clustalDir/$outputDir/;
    $bsmlFile =~ s/clustal/bsml/;

    system( "../MSF2Bsml.pl -f $clustalFile -o $bsmlFile" );
}
