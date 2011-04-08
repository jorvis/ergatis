#!/usr/bin/perl
use strict;
use warnings;

#******************************************************************************
# *merge_rarefaction_data.pl
# Author: james robert white, james.dna.white@gmail.com

# This function takes QIIME formatted mapping file and 
# a list of rarefaction curve data from mothur and generates
# csv files for input to leech. 
# Outputs are form: prefix.< mappinglevel >.csv given directly to leech
#******************************************************************************
use Getopt::Std;
use Data::Dumper;

use vars qw/$opt_f $opt_m $opt_p/;

getopts("f:m:p:");

my $usage = "Usage:  $0 \
                -f rarefaction curve file list separated by sample \
                -m mapping file from Qiime \
                -p output prefix for the csv files\
                \n";

die $usage unless defined $opt_f
              and defined $opt_m
              and defined $opt_p;

my $rarefactionfilelist = $opt_f;
my $mapfile = $opt_m;
my $prefix  = $opt_p;


# to begin, load up the mapping file levels
my %mapData = ();
loadMapData();

# now isolate the 0.03 curves generated for each sample
# SampleIDs should match perfectly
# two vectors need to be stored for each sample x,y
my %rareData = ();
my $longestseries = 0; # this will hold the length of the longest array
# how many rarefaction files are provided?
my $catstr   = `cat $rarefactionfilelist`;
my @catstr   = split "\n", $catstr;
my $numfiles = $#catstr+1;

# for each file...
for my $i (0 .. $#catstr){

  # get the filename within the directory
  my @line = split /\//, $catstr[$i];
  my $filename = $line[$#line];
  # predict SampleID
  my @b = split /\./, $filename;
  my $sampleid = join(".", @b[1..($#b-1)]);
  
  open IN, "$catstr[$i]" or die "Can't open $catstr[$i]!!\n"; 
  my $c = 0;
  while(<IN>){ # grab first two columns
    chomp($_);
    next if ($_ =~ /numsampled/);
    my @a = split "\t", $_;
    $c++;
    push @{$rareData{$sampleid}{"X"}}, $a[0];
    push @{$rareData{$sampleid}{"Y"}}, $a[1];
  }
  
  if ($c > $longestseries){
    $longestseries = $c;
  }

} # end of files

# next for each mapping level, print out
# a csv file describing the membership of 
# each sample.
open LIST, ">$prefix.csv.list" or die "Can't open $prefix.csv.list for writing!!\n";
foreach my $k (keys %mapData){
  next if ($k =~ /(Description|BarcodeSequence|LinkerPrimerSequence)/);
  printCSVTable($k);
}


#***********************************************************************
# print out comma separate variable tables for leech 
sub printCSVTable
{
  my ($antn) = @_;
  
  print LIST "$prefix.$antn.csv\n";  

  open OUT, ">$prefix.$antn.csv" or die "Cannot open $prefix.$antn.csv for writing!!\n";
  # header line
  my @orderedids = sort keys %rareData;
  print OUT "\"count\",\"$mapData{$antn}{$orderedids[0]}\"";
  for my $i (1 .. $#orderedids){
    print OUT ",\"count\",\"$mapData{$antn}{$orderedids[$i]}\"";    
  }
  print OUT "\n";

  # now for the data
  for my $j (0 .. ($longestseries-1)){
    if (defined(@{$rareData{$orderedids[0]}{"X"}}[$j])){
      print OUT "@{$rareData{$orderedids[0]}{X}}[$j]"  
    }
    print OUT ",";
    if (defined(@{$rareData{$orderedids[0]}{"Y"}}[$j])){
      print OUT "@{$rareData{$orderedids[0]}{Y}}[$j]";
    }    
  
    for my $i (1 .. $#orderedids){
      print OUT ",";
      if (defined(@{$rareData{$orderedids[$i]}{"X"}}[$j])){
        print OUT "@{$rareData{$orderedids[$i]}{X}}[$j]"
      }
      print OUT ",";
      if (defined(@{$rareData{$orderedids[$i]}{"Y"}}[$j])){
        print OUT "@{$rareData{$orderedids[$i]}{Y}}[$j]";
      }
    }
    print OUT "\n";
  } 

  close OUT;
}



sub loadMapData
{
  # this mapping file is Qiime formatted
  # e.g.  #SampleID BarcodeSequence LinkerPrimerSequence SampleName Treatment1_p Treatment2_p Descr Description
  my @maporder = ();
 
  open IN, "$mapfile" or die "Cannot open $mapfile!\n";
  while(<IN>){
    chomp($_);
    if ($_ =~ /^#/){
      my @A = split "\t", $_;
      $maporder[0] = "SampleID";
      for my $i (1 .. $#A){
        push @maporder, $A[$i];
      }
    }else{
      my @A = split "\t", $_;
      for my $i (0 .. $#A){
        $mapData{$maporder[$i]}{$A[0]} = $A[$i];
      }
    }
  }
  close IN;
#  print Dumper(\@maporder);
# print Dumper(\%mapData);
}







