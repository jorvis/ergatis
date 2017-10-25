#!/usr/bin/env perl

eval 'exec /usr/bin/env perl  -S $0 ${1+"$@"}'
    if 0; # not running under some shell
use strict;
use warnings;
use Data::Dumper;

my $wc_l = `wc -l $ARGV[0] | cut -f1 -d' '`;	# number of lines = number of OTUs present

open FILE, "$ARGV[0]"; # uclust picked otus seqs.screened_otus.txt 
my @lines = <FILE>; 
close FILE;

open OUT, ">$ARGV[2]" or die; # qiime_to_otus.txt 
open OUT2, ">$ARGV[2].list" or die; # list file
print OUT2 "$ARGV[2]"; 
close OUT2;

print OUT "$ARGV[1]\t$wc_l"; # OTU distance, numberOTUs 

my %goodseqs = ();
for my $i (0 .. $#lines){ 
  print OUT "\t"; 
  my @B = split "\t", $lines[$i]; 
  chomp($B[$#B]); 
  print OUT "$B[1]";
  $goodseqs{$B[1]} = 1;
  
  if ($#B >= 2){
    for my $j (2 .. $#B){
      $goodseqs{$B[$j]} = 1; 
      print OUT ",$B[$j]"; # printing OTUs
    } 
  }
} 
print OUT "\n";
close OUT;

open FILE, "$ARGV[3]" or die; # mothur.groups
open OUT, ">$ARGV[4]" or die; # qiime_to_otulist.groups
open OUT2, ">$ARGV[4].list" or die; # list fi;e
print OUT2 "$ARGV[4]";
close OUT2;

while(<FILE>){
  chomp($_);
  my @C = split "\t", $_;
  $C[0] =~ s/\s//g;
  if (defined($goodseqs{$C[0]})){ # keep it
    print OUT "$_\n";
  }
}
close FILE;
close OUT;



