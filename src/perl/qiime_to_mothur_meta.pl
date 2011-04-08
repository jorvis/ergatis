#!/usr/bin/perl

eval 'exec /usr/bin/perl  -S $0 ${1+"$@"}'
    if 0; # not running under some shell
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

use vars qw/$opt_f $opt_m $opt_p/;

getopts("f:m:p:");

my $usage = "Usage:  $0 \
                -f seqs.fna file from split libraries in Qiime\
                -m mapping file from Qiime \
                -p output prefix for the mothur files\
                \n";

die $usage unless defined $opt_f
              and defined $opt_m
              and defined $opt_p;

my $splitlibsfile = $opt_f;
my $mapfile = $opt_m;
my $prefix  = $opt_p;

open GROUP, ">$prefix.groups" or die;
#create list file
`echo $prefix.groups >$prefix.groups.list`;
open IN, "$splitlibsfile" or die "Can't open $splitlibsfile!!\n";
while(<IN>){
  if ($_ =~ /^>/){
    chomp($_);
    my @A = split " ", $_;
    my $name  = substr($A[0],1);
    my @namesplit = split /\_/, $name;
    my $gr = join("_", @namesplit[0..($#namesplit-1)]);
    print GROUP "$name\t$gr\n";
  }
}
close IN;
close GROUP;


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
      $pri = 1;  
    }
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
close META;

