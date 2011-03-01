#!/usr/bin/perl

eval 'exec /usr/bin/perl  -S $0 ${1+"$@"}'
    if 0; # not running under some shell
use strict;
use warnings;

#******************************************************************************
# *clovr_metagenomics_prep.pl
# Author: james robert white, james.dna.white@gmail.com

# This function takes the input parameters for the clovr metagenomics pipeline 
# and determines how to preprocess the data for the pipeline.

# A list of 1 or more fasta files is given
# Samples are described by the associated meta file. 
# These must be processed for the pipeline.
# 
# The outputs of this program are *.processed.fna and *.processed.map
#******************************************************************************
use Getopt::Std;
use File::Copy;

use vars qw/$opt_f $opt_m $opt_p/;

getopts("f:m:p:");

my $usage = "Usage:  $0 \
                -f list of fasta input files\
                -m input mapping file\
                -p output prefix for the processed files\
                \n";

die $usage unless defined $opt_f
              and defined $opt_m
              and defined $opt_p;
#
my $list    = $opt_f;
my $mapfile = $opt_m;
my $prefix  = $opt_p;
#
my $finalseqfile = "$prefix.processed.fasta"; 
my $finalmapfile = "$prefix.processed.map";
my $finalreadmapfile = "$prefix.processed.readmap.txt";
my $endlistfile  = "$prefix.processed.list";
# for now we're just copying the map file
if(!(-e $mapfile)){
  print STDERR "$mapfile does not exist! Stopping ...\n";
  exit(1);
}
copy($mapfile, $finalmapfile);

# how many fasta files are provided?
my $listlength = `wc $list`;
my @listlength = split " ", $listlength;

# create list file for download tag
open OUT, ">$endlistfile" or die "Can't open $endlistfile!!\n";
print OUT "$finalreadmapfile";
close OUT;

# we assume that each file represents a different specific sample
my $catstr = `cat $list`;
my @catstr = split "\n", $catstr;
my %filenames = ();
open SEQ, ">$finalseqfile" or die;
open READMAP, ">$finalreadmapfile" or die;
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


  # do some processing of the filename to get the prefix
  # and store the associated barcode
  my @line = split /\//, $catstr[$i];
  my $filename = $line[$#line];

  #open and process this file
  open IN, "$catstr[$i]" or die "Can't open $catstr[$i] for preprocessing!\n";
  my $seqcount = 1;
  while(<IN>){
    chomp($_);
    if ($_ =~ /^>/){ #we've reached a new sequence in this fasta file
      my @A = split " ", $_;
      my $seqname = substr($A[0],1);
      print READMAP "$filename :: $seqname :: $filename\_$seqcount\n";
      print SEQ ">$filename\_$seqcount\n";
      $seqcount++;
    }else{ 
      print SEQ "$_\n";
    }
  }
  close IN;
  # end of this file
  # move on to the next file
} 
#end of all files
close SEQ;
close READMAP;



