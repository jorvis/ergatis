#! /local/perl/bin/perl

use strict;
use warnings;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);

my %options = ();
my $results = GetOptions( \%options, 'fastaDir|f=s', 'outputDir|o=s' );

my $fastaDir = $options{'fastaDir'};
my $outputDir = $options{'outputDir'};

# remove trailing slashes from the directory names if present

$outputDir =~ s/\/+$//;
$fastaDir =~ s/\/+$//; 

die "No fasta directory specified\n" if( !$fastaDir || ! -d $fastaDir );

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

foreach my $fastaFile ( <$fastaDir/*.fasta> )
{
    my $clustalFile = $fastaFile;
    $clustalFile =~ s/$fastaDir/$outputDir/;
    $clustalFile =~ s/fasta/clustal/;

    # clustalw spits out a lot of stuff to STDOUT this redirects it to /dev/null
    system( "clustalw -output=gcg -infile=$fastaFile -outfile=$clustalFile > /dev/null" );

    unlink "$fastaDir/*.dnd";
}
