#!/usr/bin/perl
use strict;
use warnings;

#******************************************************************************
# *dataprep_dual16S.pl
# Author: james robert white, james.dna.white@gmail.com

# This function takes the input parameters for the 16S pipeline and 
# determines how to preprocess the data for the pipeline.

# If a single fasta is provided we assume it's a barcoded sample. Don't 
# change anything to it or the mapping file. 

# If multiple fasta files are given we assume the sequences are trimmed and 
# each individual fasta file is a different sample. Samples are described
# by the associated meta file. These must be processed for the pipeline.
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

my $list    = $opt_f;
my $mapfile = $opt_m;
my $prefix  = $opt_p;

my $finalseqfile = "$prefix.processed.fasta"; 
my $finalmapfile = "$prefix.processed.map";

my %data = ();
my $ARTIFICIAL_PRIMER = "CATGCTGCCTCCCGTAGGAGT";
my @barcodes = qw/AAAAAAAA AAAACGCA AAAATAGA 
                  AAACAGTA AAACGCAA AAACTTCA AAAGCCGA AAAGGTTA AAATAGAA AAATGACA AAATTGGA AACACATA 
                  AACAGTAA AACCACCA AACCCTGA AACCTCTA AACGCAAA AACGGGCA AACTAAGA AACTCGTA AACTTCAA 
                  AATCGGTA AATGACAA AATGCTCA AATGTCGA AATTATTA AATTGGAA ACAAAACA ACAACGGA ACAATATA ACACATAA ACACGCCA ACACTTGA ACAGCCTA ACAGTAAA 
                  ACATAGCA ACATGAGA ACATTGTA ACCACCAA ACCAGTCA ACCCACGA ACCCCTTA ACCCTGAA ACCGCACA ACCGGGGA ACCTAATA ACCTCTAA ACCTTCCA ACGAATGA 
                  ACGAGCTA ACGCAAAA ACGCCGCA ACGCTAGA ACGGAGTA ACGGGCAA ACGGTTCA ACGTCCGA ACGTGTTA ACTAAGAA ACTAGACA ACTATGGA ACTCCATA ACTCGTAA 
                  AGTTGGGA ATAAAATA ATAACTAA ATAATCCA ATGAGGCA ATGCAAGA ATGCCGTA ATGCTCAA ATGGATCA ATGGGCGA ATGGTTTA ATGTCGAA ATGTTACA ATTAAGGA 
                  ATTAGATA ATTATTAA ATTCCCCA ATTCGTGA ATTGACTA ATTGGAAA ATTGTGCA ATTTCAGA ATTTGGTA CAAAACAA CAAACTCA CAAATCGA CAACATTA CAACGGAA 
                  CAAGAACA CAAGCGGA CAAGTATA CAATATAA CAATGCCA CAATTTGA CACACCTA CACATAAA CACCAGCA CACCGAGA CACCTGTA CACGCCAA CACGGTCA CACTACGA 
                  CACTCTTA CACTTGAA CAGACACA CAGAGGGA CAGCAATA CAGCCTAA CAGCTCCA CAGGATGA CAGGGCTA CAGTAAAA CAGTCGCA CAGTTAGA CATAAGTA CATAGCAA 
                  CATATTCA CATCCCGA CATCGTTA CATGAGAA CATGGACA CATGTGGA CATTCATA CATTGTAA CCAAACCA CCAACTGA CCAATCTA CCACCAAA CCACGGCA CCAGAAGA 
                  CCAGCGTA CCAGTCAA CCATATCA CCATGCGA CCATTTTA/;

# how many fasta files are provided?
my $listlength = `wc $list`;
my @listlength = split " ", $listlength;
if ($listlength[0] > $#barcodes-1){
  print STDERR "We have not implemented enough artificial barcodes to properly handle your multple fasta file input to the pipeline.\nSorry!\n"; 
  exit(1);
}
if ($listlength[0] == 1){
  my $file = `cat $list`;
  chomp($file);
  if((-e $file) and (-e $mapfile)){
    copy($file, $finalseqfile);
    copy($mapfile, $finalmapfile);
    print STDOUT "One fasta file detected. We assume this file is barcoded to determine samples and also that the mapping file provided in formatted for Qiime.\n"; 
    exit(0);
  }elsif(!(-e $file)){
    print STDERR "$file does not exist! Stopping ...\n";
    exit(1);
  }elsif(!(-e $mapfile)){
    print STDERR "$mapfile does not exist! Stopping ...\n";
    exit(1);
  }
}

# otherwise we've got multiple fasta files
# we assume that these files are trimmed of barcodes/primers and that each file
# represents a different specific sample.
my $catstr = `cat $list`;
my @catstr = split "\n", $catstr;
open SEQ, ">$finalseqfile" or die;
for my $i (0 .. ($listlength[0]-1)){
  my $bc = $barcodes[$i]; # this is the associated barcode for the sample

  # do some processing of the filename to get the prefix
  # and store the associated barcode
  my @line = split /\//, $catstr[$i];
  $data{$line[$#line]} = $bc;
  my @linesplit = split /\./, $line[$#line];
  my $fileprefix = join(".", @linesplit[0..($#linesplit-1)]);

  if(!(-e $catstr[$i])){
    print STDERR "$catstr[$i] does not exist! Check your file names. Stopping ...\n";
    exit(1);
  }

  #open and process this file
  open IN, "$catstr[$i]" or die "Can't open $catstr[$i] for preprocessing!\n";
  my $seq = "";
  my $seqname = "";
  my $seqcount = 1;
  my $fL = 70;
  while(<IN>){
    chomp($_);
    if ($_ =~ /^>/){ #we've reached a new sequence in this fasta file
      my @A = split " ", $_;
      if ($seq ne ""){
        print SEQ ">$fileprefix\_$seqcount\n";
        my $tmp  = $bc;
        $tmp    .= $ARTIFICIAL_PRIMER;
        $tmp    .= $seq;
 
        for (my $i = 0; $i<length($tmp); $i+=$fL){
          my $substr = substr($tmp, $i, $fL);
          print SEQ "$substr\n";
        }
        $seqcount++;
      }
      $seqname = substr($A[0],1);
      $seq = "";
    }else{ 
      $seq .= "$_";
    }
  }
  close IN;

  print SEQ ">$fileprefix\_$seqcount\n";    
  my $tmp = $bc;
  $tmp .= $ARTIFICIAL_PRIMER;
  $tmp .= $seq;

  for (my $i = 0; $i<length($tmp); $i+=$fL){
    my $substr = substr($tmp, $i, $fL);
    print SEQ "$substr\n";
  }
  # end of this file

}
close SEQ;

# now create the corresponding mapping file
if(!(-e $mapfile)){
  print STDERR "$mapfile does not exist! Stopping ...\n";
  exit(1);
}

open MAP, ">$finalmapfile" or die;
open IMAP, "$mapfile" or die "Can't open $mapfile!!\n";
while(<IMAP>){
  chomp($_);
  my @A = split " ", $_;
  if ($_ =~ /^#/){
    print MAP "#SampleID\tBarcodeSequence\tLinkerPrimerSequence";
    for my $j (1 .. $#A){
      print MAP "\t$A[$j]";
    }
    print MAP "\n";
  }else{
    print MAP "$A[0]\t$data{$A[0]}\t$ARTIFICIAL_PRIMER";
    for my $j (1 .. $#A){
      print MAP "\t$A[$j]";
    }
    print MAP "\n";
  }
}
close IMAP;
close MAP;



