#!/usr/bin/perl
use strict;
use warnings;

#******************************************************************************
# *qiime_to_mothur_meta.pl
# Author: james robert white, james.dna.white@gmail.com

# This function takes QIIME formatted mapping file and 
# converts it to an oligos file and a metadata file
# for the mothur component of the dual 16S pipeline
#******************************************************************************
use Getopt::Std;
use warnings;

use vars qw/$opt_m $opt_p/;

getopts("m:p:");

my $usage = "Usage:  $0 \
                -m qiime formatted mapping file\
                -p output prefix for the oligos and meta files\
                \n";

die $usage unless defined $opt_m
              and defined $opt_p;

my $mapfile = $opt_m;
my $prefix  = $opt_p;

open OLI, ">$prefix.oligos" or die;
open META, ">$prefix.meta" or die;
open IN, "$mapfile" or die;
my $pri = 0;
my @pairwise_indices = ();
while(<IN>){
  chomp($_);
  if ($_ =~ /^#SampleID/){
    my @A = split "\t", $_;
    for my $i (0 .. $#A){
      if (substr($A[$i], -2) eq "_p"){
        $pairwise_indices[$i] = 1;    
      }else{
        $pairwise_indices[$i] = 0;
      } 
    }   
    next;
  }else{
    my @A = split "\t", $_;
    if ($pri == 0){
      print OLI "forward\t$A[2]\n";
      $pri = 1;  
    }
    print OLI  "barcode\t$A[1]\t$A[0]\n";
    print META "barcode, $A[1], $A[0]";
    if ($#A > 2){
      for my $i (3 .. $#A){
        if ($pairwise_indices[$i] ==1){
          print META ", $A[$i]";
        }  
      }
    }
    print META "\n";
  }  
}
close IN;


