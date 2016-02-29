#!/usr/bin/perl
use strict;
use warnings;

#******************************************************************************
# *clovr_map_MiSeq_check.pl
# Author: Shaun Adkins - sadkins@som.umaryland.edu
#
# This function takes the input parameters for the clovr 16S pipeline 
# and determines if they are consistent enough for a realistic run.

# A list of 2 fastq files must be given
# Samples are described by the associated mapping file.
# 
# There are no file outputs in the scripts. It either passes or fails.
#******************************************************************************
use Getopt::Std;
use vars qw/$opt_f $opt_m $opt_b/;

getopts("f:m:b:");

my $usage = "Usage:  $0 \
                -f list of fastq input files\
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

# how many fastq files are provided?
my $listlength = `wc $list`;
my @listlength = split " ", $listlength;

if ($listlength[0] == 0){
  print STDERR "Empty fastq filelist detected. Need to tag sequence data. Stopping ...\n";
  exit(1);
}

# if it's qiime then there should only be one fastq file
if ($listlength[0] != 2){
  print STDERR "For FastQ input, exactly 2 FastQ sequence files must be provided. This is not currently compatible with the pipeline. Stopping ...\n";
  exit(1);
}

exit(0);

