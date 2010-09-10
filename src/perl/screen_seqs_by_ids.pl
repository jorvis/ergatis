#!/usr/bin/perl

eval 'exec /usr/bin/perl  -S $0 ${1+"$@"}'
    if 0; # not running under some shell
use strict;
use warnings;

#******************************************************************************
# *screen_seqs_by_ids.pl
# Author: james robert white, james.dna.white@gmail.com

# This function takes a list of fasta files and a list of sequence ids as input.
# Depending on the parameters set by the user, this function will either output
# a fasta file containing all seqs with the ids in the list, or the set that
# are not in the list.
#******************************************************************************
use Getopt::Std;
use Data::Dumper;
use warnings;

use vars qw/$opt_f $opt_s $opt_o $opt_r/;

getopts("f:s:o:r");

my $usage = "Usage:  $0 \
                -f <list of fasta files>\
                -s <list of sequence ids>\
                -o <output file name>\
                -r set this flag to output all seqs
                   NOT in the list 
                \n";

die $usage unless defined $opt_f
              and defined $opt_s
              and defined $opt_o;

# Vars:
my %seqids     = ();
my $outputfile = $opt_o;
my $fastalist  = $opt_f;
my $seq_list   = $opt_s;
my @fastafiles = ();

# first load up the seq_ids
loadIds();

# now open up the fasta list and process the files
open OUT, ">$outputfile" or die "Cannot open $outputfile!!\n";
if (defined($fastalist)) {
  open (FLIST, $fastalist) or die ("Could not open FASTA list $fastalist: $!");
  while (my $file = <FLIST>) {
    chomp($file);
    push (@fastafiles, $file);
  }
  close (FLIST);
}


foreach my $file (@fastafiles){  
  open IN, $file or die "Could not open FASTA file $file!!\n";
  my $ck = 0;
  while(<IN>){
    chomp($_);
    if ($_ =~ /^>/){
      my @A = split " ", $_; 
      my $name = substr($A[0],1);
             
      if (!$opt_r){
        if (defined($seqids{$name})){
          $ck = 1;
        }else{
          $ck = 0;
        } 
      }else{
        if (defined($seqids{$name})){
          $ck = 0;
        }else{
          $ck = 1;
        }
      }
       
    }

    print OUT "$_\n" if ($ck == 1);

  }
  close IN;
}

close OUT;








sub loadIds
{
  open IN, "$seq_list" or die "Cannot open $seq_list!!\n";
  while(<IN>){
    chomp($_);
    $seqids{$_} = 1; 
  }
  close IN;

}
 





