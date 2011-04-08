#!/usr/bin/perl
use strict;
use warnings;

#******************************************************************************
# *merge_mothur_chimera_ids.pl
# Author: james robert white, james.dna.white@gmail.com

# This function takes a list of files with chimera sequence ids and a list of 
# files of clusters from mothur_unique_seqs as input.

# The output is a single concatenated list of ids that are all chimeric 
# sequences.
#******************************************************************************
use Getopt::Std;
use Data::Dumper;

use vars qw/$opt_n $opt_c $opt_o/;

getopts("n:c:o:");

my $usage = "Usage:  $0 \
                -n <list of name cluster files>\
                -c <list of files with representative chimera sequence ids>\
                -o <output file name>\
                \n";

die $usage unless defined $opt_n
              and defined $opt_c
              and defined $opt_o;

my $uniqueRepFileList = $opt_n;
my $chimeraFileList   = $opt_c;
my $outfile           = $opt_o;

my @uniqueclusterfiles = ();
my %clusters = ();

# load up list of unique read clusters
open (FLIST, $uniqueRepFileList) or die ("Could not open Unique read file list $uniqueRepFileList: $!");
while (my $file = <FLIST>) {
  chomp($file);
  push (@uniqueclusterfiles, $file);
}
close (FLIST);

# create a hash of all clusters
foreach my $file (@uniqueclusterfiles){
  open IN, $file or die "Could not open file $file!!\n";
  while(<IN>){
    chomp($_);
    my @A = split "\t", $_;
    my @B = split ",", $A[1];
    for my $i (0 .. $#B){
      push @{$clusters{$B[0]}}, $B[$i];
    }# this cluster includes the representative in the hash value
  }
}

# now go through the list of chimera accnos and 
# print out all chimera ids using the clusters hash
my @chifiles = ();
open (FLIST, $chimeraFileList) or die ("Could not open Unique read file list $chimeraFileList: $!");
while (my $file = <FLIST>) {
  chomp($file);
  push (@chifiles, $file);
}
close (FLIST);


open OUT, ">$outfile" or die;

foreach my $file (@chifiles){
  open IN, $file or die "Could not open file $file!!\n";
  while(<IN>){
    chomp($_);
    if (!defined($clusters{$_})){
      die "Cannot locate representative sequence: $_!!\n"; 
    }else{ #print that cluster out b/c theyre all chimeras
      foreach my $v (@{$clusters{$_}}){
        print OUT "$v\n"; 
      }
    } 
  }
  close IN;
}
close OUT;


