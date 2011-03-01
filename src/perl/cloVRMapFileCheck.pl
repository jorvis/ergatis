#!/usr/bin/perl
use strict;
use warnings;

# author: james.dna.white@gmail.com
# Created: 03-01-2011
# This script checks correctness of a submitted mapping file
# in either Qiime of CloVR map format.

# Will error with description if format is wrong.

# First check the mapping file type
my $file     = $ARGV[0]; # the mapping file is provided as the first argument on the command line
my $filetype = "";       # this will end up being qiime or clovr
open IN, "$file" or die "Can't open $file for reading!!\n";
my $L = 0; # Line number count
while(<IN>){
  chomp($_);
  my @A  = split "\t", $_;
  if ($L == 0){  
    if ($A[0] eq "\#File"){        # then it's a clovr format mapping
      $filetype = "clovr";
    }elsif($A[0] eq "\#SampleID"){ # then it's a qiime format mapping
      $filetype = "qiime";
    }else{
      die "Cannot determine the type of mapping file! Header line must begin with #File or #SampleID\n";
    }
  }    
  
  $L++;
  my @B = split "", $_;

  # No empty lines!
  if (!defined($B[0])){
    die "Error: Please remove any empty lines from your mapping file.\n";
  }
  
  # No spaces! Tabs only!
  for my $i (0 .. $#B){
    if ($B[$i] eq " "){
      die "Whitespace found on line $L of mapping file. Will likely cause an error. Please remove all non-tab whitespace before proceeding.\n";
    }
  }  
}
close IN;

if ($filetype eq "qiime"){
  checkQiimeFormat();  
}elsif($filetype eq "clovr"){
  checkCloVRFormat();
}


#*******************************************************************
sub checkQiimeFormat
{
  my $l = 0;
  print "Checking Qiime-formatted mapping file...\n";
  open IN, "$file" or die "Can't open $file for reading!!\n";
  my %samplenames  = ();
  my $headerlength = ();
  while(<IN>){
    chomp($_);
    my @A = split "\t", $_;
     
    if ($l == 0){ # then this is the header line, check reqs   
      if ($A[0] ne "\#SampleID"){
        die "Error: Your Qiime-formatted mapping file must begin with \"#SampleID\". Remember to remove all white spaces.\n";
      }
      if ($A[1] ne "BarcodeSequence"){
        die "Error: The second column header must be \"BarcodeSequence\" in a Qiime-formatted mapping file. Remember to remove all white spaces.\n";
      }
      if ($A[2] ne "LinkerPrimerSequence"){
        die "Error: The third column header must be \"LinkerPrimerSequence\" in a Qiime-formatted mapping file. Remember to remove all white spaces.\n";
      }
      if ($A[$#A] ne "Description"){
        die "Error: The final column header must be \"Description\" in a Qiime-formatted mapping file. Remember to remove all white spaces.\n";
      }
         
      $headerlength = $#A+1;
      if ($#A > 3){
        for my $i (3 .. ($#A-1)){
          if ($A[$i] !~ /[A-Za-z0-9\_\.\-\+\%]+/ or $A[$i] =~ /[\?\!\@\$\^\&\*\(\)\'\"\;\:\<\>]+/){
            die "Error: A column header must be a combination of only alphanumeric, underscore (_), period (.), minus sign (-), plus sign (+) and/or percentage (%) characters.";
          }
        }
      }
      $l++;
      next;

    }elsif($l>0){ # then we're looking below the header line, so check data types
    
      if (defined($samplenames{$A[0]})){
        die "Error: SampleID $A[0] is defined multiple times in the mapping file!\n";
      }else{
        $samplenames{$A[0]} = 1;
        if ($A[0] !~ /[A-Za-z0-9\.]+/ or $A[0] =~ /[\?\!\@\$\^\&\*\(\)\'\"\;\:\<\>]+/){
          die "Error: $A[0] is an invalid sampleid for a qiime-formatted mapping file. Please use alphanumerics and period characters.\n"; 
        }
      } 

      if ($headerlength != ($#A+1)){
        die "Error: Row $l of the mapping file has a different number of elements than the header line.\n";
      }else{
        for my $i (0 .. $#A){
          my $j = $i+1;
          if (!defined($A[$i])){
            die "Error: In row $l of mapping file, the entry in column $j is undefined.\n"; 
          }
          if ($A[$i] !~ /[A-Za-z0-9\_\.\-\+\%]+/ or $A[$i] =~ /[\?\!\@\$\^\&\*\(\)\'\"\;\:\<\>]+/){
            die "Error: In row $l of mapping file, the entry in column $j is not a combination of only alphanumeric, underscore (_), period (.), minus sign (-), plus sign (+) and/or percentage (%) characters.";
          } 
        } 
      }
      $l++;
      next;
    }
  }
  close IN;
  print "Your Qiime-formatted mapping file looks good!\n";

}


sub checkCloVRFormat
{
  
  my $l = 0;
  print "Checking CloVR-formatted mapping file...\n";
  open IN, "$file" or die "Can't open $file for reading!!\n";
  my %samplenames  = ();
  my $headerlength = ();
  while(<IN>){
    chomp($_);
    my @A = split "\t", $_;

    if ($l == 0){ # then this is the header line, check reqs   
      if ($A[0] ne "\#File"){
        die "Error: Your CloVR-formatted mapping file must begin with \"#File\". Remember to remove all white spaces.\n";
      }
      if ($A[1] ne "SampleName"){
        die "Error: The second column header must be \"SampleName\" in a CloVR-formatted mapping file. Remember to remove all white spaces.\n";
      }

      $headerlength = $#A+1;
      for my $i (2 .. $#A){
        if ($A[$i] !~ /[A-Za-z0-9\_\.\-\+]+/ or $A[$i] =~ /[\?\!\@\$\^\&\*\(\)\'\"\;\:\<\>]+/){
          die "Error: A column header must be a combination of only alphanumeric, underscore (_), period (.), minus sign (-), plus sign (+) characters.";
        }
      }
      $l++;
      next;

    }elsif($l>0){ # now we are below the header line 

      if (defined($samplenames{$A[0]})){
        die "Error: SampleID $A[0] is defined multiple times in the mapping file!\n";
      }else{
        $samplenames{$A[0]} = 1;
        if ($A[0] !~ /[A-Za-z0-9\.\_]+/ or $A[0] =~ /[\?\!\@\$\^\&\*\(\)\'\"\;\:\<\>\%]+/ ){
          die "Error: $A[0] is an invalid sampleid for a qiime-formatted mapping file. Please use alphanumeric, underscore, and period characters.\n";
        }
      }

      if ($headerlength != ($#A+1)){
        die "Error: Row $l of the mapping file has a different number of elements than the header line.\n";
      }else{
        for my $i (0 .. $#A){
          my $j = $i+1;
          if (!defined($A[$i])){
            die "Error: In row $l of mapping file, the entry in column $j is undefined.\n";
          }
          if ($A[$i] !~ /[A-Za-z0-9\_\.\-\+\%]+/ or $A[$i] =~ /[\?\!\@\$\^\&\*\(\)\'\"\;\:\<\>]+/){
            die "Error: In row $l of mapping file, the entry in column $j is not a combination of only alphanumeric, underscore (_), period (.), minus sign (-), plus sign (+) and/or percentage (%) characters.";
          }
        }
      }
      $l++;
      next;
    }
  }
  close IN;
    print "Your CloVR-formatted mapping file looks good!\n"; 
}

