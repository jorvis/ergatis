#! /usr/local/bin/perl 

use strict;
use warnings;
use Data::Dumper;

open( INFILE, "$ARGV[0]" ) or die "Could not open COGS data\n";

my $cogs = {};
my $genomeIds = {};

while( my $line = <INFILE> )
{
    if( $line =~ /^COG = ([\d]*), size ([\d]*)/ )
    {
	my $proteins = [];
	my $cogId = $1;
	my $size = $2;
	
	for( my $i=0; $i<$size; $i++ )
	{
	    my $protLine = <INFILE>;
	    $protLine =~ /\t([\S]*)\n/;

	    push( @{$proteins}, $1 );

	    my @id = split /\./, $1;
	    $genomeIds->{$id[0]} = 1;
	}

	$cogs->{$cogId}->{'size'} = $size;
	$cogs->{$cogId}->{'proteins'} = $proteins;
    }
    else
    {
	print "BAD Line: $line\n";
	exit();
    }
}

# Find the total number of Clusters

my $cogsCount = keys(%{$cogs});

print "Total number of clusters in Analysis: $cogsCount\n";

# Find the number of clusters with two members

my $pairCount = 0;
my $ana1AOA1Count = 0;
my $ana1AFU1Count = 0;
my $aoa1AFU1Count = 0;

my $pairCounts = {};

foreach my $key (keys(%{$cogs}))
{
    if( $cogs->{$key}->{'size'} == 2 )
    {
	$pairCount++;

	my @ids = keys( %{$genomeIds} );

	for( my $i=0; $i<@ids; $i++ )
	{
	    for( my $j=($i+1); $j<@ids; $j++ )
	    {
		if( (($cogs->{$key}->{'proteins'}->[0] =~ /^$ids[$i]/) || 
		     ($cogs->{$key}->{'proteins'}->[0] =~ /^$ids[$j]/)) &&
		    (($cogs->{$key}->{'proteins'}->[1] =~ /^$ids[$i]/) ||
		     ($cogs->{$key}->{'proteins'}->[1] =~ /^$ids[$j]/)) )
		{
		    $pairCounts->{"$ids[$i] - $ids[$j]"}++;
		}
	    }
	}
    }
}

print "Number of Clusters with 2 members: $pairCount\n";

foreach my $key (keys(%{$pairCounts}))
{
    print "$key : $pairCounts->{$key}\n";
}

# Find the number of clusters with three members

my $tripleCount = 0;

foreach my $key (keys(%{$cogs}))
{
    if( $cogs->{$key}->{'size'} == 3 )
    {
	$tripleCount++;
    }
}

print "Number of Clusters with 3 members: $tripleCount\n";

# Find the number of clusters with four members 

my $quadCount = 0;

foreach my $key (keys(%{$cogs}))
{
    if( $cogs->{$key}->{'size'} == 4 )
    {
	$quadCount++;
    }
}

print "Number of Clusters with 4 members: $quadCount\n";

# Find the number of clusters with five members 

my $pentCount = 0;

foreach my $key (keys(%{$cogs}))
{
    if( $cogs->{$key}->{'size'} == 5 )
    {
	$pentCount++;
    }
}

print "Number of Clusters with 5 members: $pentCount\n";

# Find the number of clusters with more than five members 

my $bigCount = 0;

foreach my $key (keys(%{$cogs}))
{
    if( $cogs->{$key}->{'size'} > 5 )
    {
	$bigCount++;
    }
}

print "Number of Clusters with more than 5 members: $bigCount\n";

