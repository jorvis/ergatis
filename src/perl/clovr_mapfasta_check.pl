#!/usr/bin/perl
use strict;
use warnings;

#******************************************************************************
# *clovr_mapfasta_check.pl
# Author: james robert white, james.dna.white@gmail.com

# This function takes the input parameters for the clovr metagenomics pipeline
# or the clovr 16S pipeline and determines if they are consistent enough for 
# a realistic run.

# A list of 1 or more fasta files is given
# Samples are described by the associated mapping file. 
# 
# There are no file outputs in the scripts. It either passes or fails.
#******************************************************************************
use Getopt::Std;
use vars qw/$opt_f $opt_m/;

getopts("f:m:");

my $usage = "Usage:  $0 \
                -f list of fasta input files\
                -m input mapping file\
                \n";

die $usage unless defined $opt_f
              and defined $opt_m;

my $list    = $opt_f;
my $mapfile = $opt_m;
##

# does the mapping file exist?
if(!(-e $mapfile)){
  print STDERR "$mapfile does not exist! Stopping ...\n";
  exit(1);
}

my %mappingfilenames = ();
my $maptype = ""; 
open IN, "$mapfile" or die;
my $ck = 0;
while(<IN>){
  chomp($_);
  my @A = split "\t", $_;
  if ($ck == 0){
    if ($A[0] eq "#File"){
      $maptype = "clovr";
    }else{
      $maptype = "qiime";
      last;
    }
    $ck = 1;
    next;
  }else{
    $mappingfilenames{$A[0]} = 1;
  }
}
close IN;
#

# how many fasta files are provided?
my $listlength = `wc $list`;
my @listlength = split " ", $listlength;

if ($listlength[0] == 0){
  print STDERR "Empty fasta filelist detected. Need to tag sequence data. Stopping ...\n";
  exit(1);
}

# if it's qiime then there should only be one fasta file
if ($maptype eq "qiime" and $listlength[0] > 1){
  print STDERR "Mapping file is Qiime-based but there is more than one fasta file. We assume for a Qiime-based mapping file that there is one fasta file and that it has multiplex sample-specific barcodes. This is not currently compatible with the pipeline. Stopping ...\n";
  exit(1);
}

if ($maptype eq "qiime"){
  return;
} # everything below this is for clovr-formatted mapping files

# we assume that each file represents a different specific sample
my $catstr = `cat $list`;
my @catstr = split "\n", $catstr;
my %filenames = ();
for my $i (0 .. ($listlength[0]-1)){
  if(!(-e $catstr[$i])){
    print STDERR "$catstr[$i] does not exist! Check your file names. Stopping ...\n";
    exit(1);
  }
  if (!defined($filenames{$catstr[$i]})){
    $filenames{$catstr[$i]} = 1;
  }else{
    print STDERR "$catstr[$i] is being used more than once! Check your file names. Stopping ...\n";
    exit(1);
  }
} 
# end of all fasta files


# now check the consistency of the map to the fasta list
my %shortflnames = ();
foreach my $fl (keys %filenames){
  my @A = split /\//, $fl;
  $shortflnames{$A[$#A]} = 1;
  if (!defined($mappingfilenames{$A[$#A]})){
    print STDERR "Warning the fasta filelist has the file $A[$#A], but the mapping file does not. Check your filenames. Stopping ...\n";
    exit(1);
  }
}

foreach my $ml (keys %mappingfilenames){
  if (!defined($shortflnames{$ml})){
    print STDERR "Warning the mapping file has $ml, but the fasta filelist does not. Check your filenames. Stopping ...\n";
    exit(1);
  }
}


