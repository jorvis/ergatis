#!/usr/bin/perl
use strict;
use warnings;
use Data::Dumper;

#*********************************************************************
#  Create_Qiime_tables.pl*
#  * This program takes taxonomic tables created in the initial qiime 
#  pipeline and reformats them for comparative analysis with skiff and 
#  metastats.

#  The output is *.tsv and *.2tsv

#  Author: james robert white, whitej@umd.edu  
#  Last modified: August 11, 2010 
#*********************************************************************
use Getopt::Std;
use warnings;

use vars qw/$opt_t $opt_m $opt_p/;

getopts("t:m:p:");

my $usage = "Usage:  $0 \
                -t otu summary table from qiime pipeline\
                -m metadata file describing associations between samples\
                -p output prefix [<prefix>.tsv]
                \n";

die $usage unless defined $opt_m
              and defined $opt_p
              and defined $opt_t;

my $prefix = $opt_p;

#*********************************************************************
# GLOBALS
#*********************************************************************
my @samples        = ();
my %activesamples  = ();
my @orderedsamples = ();
my %groups         = ();
my %seqmap         = ();
my %data           = ();
my %activefeatures = ();
my %groupnames = ();
#*********************************************************************

#*********************************************************************
# Begin parsing files
#*********************************************************************
my $numfeaturetypes = 0;
open IN, "$opt_m" or die "Can't open Meta file for reading!!\n";
while(<IN>){
  chomp($_);
  next if ($_ eq "");
  my @A = split /\,/, $_;
  $A[2] =~s/^\s+|\s+$//g;
  $numfeaturetypes = $#A-2;
  if (!defined($A[3])){ # then there are no classes for the samples
    push @{$groups{"NoClassDef"}}, $A[2];    
  }else{
    for my $i (3 .. $#A){
      $A[$i] =~s/^\s+|\s+$//g;
      push @{$groups{$i-2}{$A[$i]}}, $A[2];  
      $groupnames{$i-2}{$A[$i]} = 1; 
    }
  }
}
close IN;

# fill up the data hash, active features, and active samples
open IN, "$opt_t" or die;
while(<IN>){
  chomp($_);
  next if ($_ =~ /^#/);
  my @A = split "\t", $_;
  if ($_ =~ /Taxon/){
    for my $i (1 ..$#A){
      $activesamples{$A[$i]} = 1;
      push @samples, $A[$i]; 
    }                
  }else{
    my @B = split /\;/, $A[0];
    my $f = $B[$#B];      
    $activefeatures{$f} = 1;
    for my $i (1 ..$#A){
      $data{$samples[$i-1]}{$f} = substr($A[$i], 0, length($A[$i])-2);  
    }        
  }
}
close IN;

if ($numfeaturetypes > 0){
foreach my $type (keys %{$groups{1}}){
  foreach my $s (@{$groups{1}{$type}}){
    if (defined($activesamples{$s})){
      push @orderedsamples, $s;
    }
  }
}
}

if (!defined($orderedsamples[0])){ 
  @orderedsamples = keys %activesamples;
}

printCounts();

for my $f (1 .. $numfeaturetypes){
  my @gs = sort keys %{$groupnames{$f}};
  for my $i (0 .. $#gs){
    for my $j (0 .. $#gs){
      next if ($j >= $i);
      printPairedGroups($f, $gs[$i], $gs[$j]);
    }
  }
}  

#*********************************************************************
# Subroutines
#*********************************************************************

sub printCounts
{
  open OUT, ">$prefix.tsv" or die "Can't open $prefix.tsv for writing!\n"; 
  
  for my $s (@orderedsamples){
    print OUT "\t$s";
  } 
   print OUT "\n";

  foreach my $f (sort keys %activefeatures){
    print OUT "$f";
    for my $s (@orderedsamples){
      print OUT "\t$data{$s}{$f}";
    } 
    print OUT "\n";
  } 
  close OUT;
}


sub printPairedGroups
{
  my ($feature, $g1, $g2) = @_;

  my @pairedorderedsamples = ();
  my $g1count = 0;
  foreach my $s (@{$groups{$feature}{$g1}}){
     if (defined($activesamples{$s})){
       push @pairedorderedsamples, $s;
       $g1count++; 
     }
  }
  my $g2count = 0;
  foreach my $s (@{$groups{$feature}{$g2}}){
     if (defined($activesamples{$s})){
     push @pairedorderedsamples, $s;
     $g2count++;
    }
  }
   
  return if ($g1count <= 0 or $g2count <=0);
  return if (($g1count == 1 and $g2count != 1) or ($g1count != 1 and $g2count == 1));  
 
  open OUT, ">$prefix.$g1\_vs_$g2.$g1count-$g2count.2tsv" or die "Can't open $prefix.$g1\_vs_$g2.$g1count--$g2count.2tsv for writing!\n";
  for my $s (0 .. $#pairedorderedsamples){
    print OUT "\t$pairedorderedsamples[$s]";
  }
  print OUT "\n";

  foreach my $f (sort keys %activefeatures){
    print OUT "$f";
    for my $s (0 .. $#pairedorderedsamples){
      print OUT "\t$data{$pairedorderedsamples[$s]}{$f}";
    }
    print OUT "\n";
  }
  close OUT;
}

