#!/usr/bin/perl
use strict;
use warnings;
use Data::Dumper;

#*********************************************************************
#  mothur_unique_seq_limit.pl*
#  Checks to see if there are too many unique sequences in the mothur
#  pipeline.
#  To mitigate pipeline workflow complaints, if there are too many 
#  sequences this program will replace the fasta set with a faux set
#  very small, everything will run to completion, but not correctly.

#  Author: james robert white, whitej@umd.edu  
#  Last modified: August 10, 2010 
#*********************************************************************
use Getopt::Std;
use File::Copy;
use warnings;

use vars qw/$opt_f $opt_n $opt_e $opt_g $opt_m/;

getopts("f:n:e:g:m:");

my $usage = "Usage:  $0 \
                -f fasta file
                -g names file
                -n sequence no. limit (e.g. 50000)
                -e faux file path (e.g. /opt/opt-packages/bioinf-v1r4b1/mothur/emptyalign.fasta)
		-m faux names path (e.g. /opt/opt-packages/bioinf-v1r4b1/mothur/empty.names)
                \n";

die $usage unless defined $opt_f and
                  defined $opt_n and
		  defined $opt_e;

my $str = `grep '>' -c $opt_f`;
chomp($str);

my $empty = $opt_e;
my $emptynames = $opt_m;
my @A = split " ", $str;
my $numseqs = $A[0];
my $N = $opt_n;

if ($numseqs > $N){
  print STDERR "\n\n**Error: there are too many unique sequences to continue the mothur pipeline.\nThe remainder of the mothur pipeline will finish using a small faux fasta file.\nSorry!\n\n"; 
  copy($empty, $opt_f); 
  copy($emptynames, $opt_g); 
}



