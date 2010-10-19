#!/usr/bin/perl

eval 'exec /usr/bin/perl  -S $0 ${1+"$@"}'
    if 0; # not running under some shell
use strict;
use warnings;

#******************************************************************************
# *uclust_replicate_removal.pl
# Author: james robert white, james.dna.white@gmail.com

# Filters a clstr file from UCLUST module 
# Removes artificial pyrosequencing replicates
#******************************************************************************
use Getopt::Std;
use File::Copy;

use vars qw/$opt_f $opt_l $opt_k $opt_p/;

getopts("f:l:k:p:");

my $usage = "Usage:  $0 \
                -f clstr file from uclust\
                -l length difference maximum\
                -k prefix length match requirement\
                -p output prefix for the processed files\
                \n";

die $usage unless defined $opt_f
              and defined $opt_l
              and defined $opt_k
              and defined $opt_p;
#
my $input    = $opt_f;
my $length  = $opt_l;
my $k       = $opt_k;
my $prefix  = $opt_p;
#
my $finalclstrfile = "$prefix.clstr"; 
my $finalrepfile   = "$prefix.replicates";

my %prefixes = ();

# first grab all associated prefixes of length $k
my @A = split /\./, $input; 
my $tmpstr = join(".", @A[0..($#A-1)]);
my $sorted_fasta = "$tmpstr.sorted";
open IN, "$sorted_fasta" or die "Can't open $sorted_fasta for prefix counting!\n";
my $cseq = "";
while(<IN>){
  chomp($_);
  if ($_ =~ /^>/){
    $cseq = $_;
  }else{
    $prefixes{$cseq} = substr($_,0,$k);
  }
}
close IN;


# now filter the clstr file for prefixes matches
# within a specific length

open IN, "$input" or die "Cannot open $input for reading!\n";
open OUT, ">$finalclstrfile" or die "Cannot open $finalclstrfile for writing!\n";
open REP, ">$finalrepfile" or die "Cannot open $finalrepfile for writing!\n";

my $repseq    = "";
my $replength = "";
while(<IN>){
  chomp($_);
  if ($_ =~ /^>/){
    print OUT "$_\n";
    next;   
  }else{
    my @B = split "nt, ", $_;
    my @C = split /\.\.\./, $B[1];
    my @D = split " ", $B[0];

    if ($_ =~ /^0/){
      $repseq = $C[0];
      $replength = $D[1];
      print OUT "$_\n";    
      next;
    }elsif(($prefixes{$repseq} eq $prefixes{$C[0]}) and (($replength - $D[1])  < $length)){
      print REP "$C[0] is a replicate of $repseq\n";
      next;
    }else{
      print OUT "$_\n";
      next;  
    } 
  }
}
close IN;
close REP;
close OUT;

