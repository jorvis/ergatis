#!/usr/bin/perl
use strict;
use warnings;

#******************************************************************************
# *create_feature_tables.pl
# Author: james robert white, james.dna.white@gmail.com

# This function takesa list blast outputs (-m 8) and some metadata and creates
# a set of files for postprocessing in the clovr metagenomics pipeline.

# A list of 1 or more blast output files is given
# Samples are described by the associated mapping file. Hits are described by
# the annotation file provided as input.
# 
# The outputs of this program are tab-delimited tables of hits
#******************************************************************************
use Getopt::Std;
use Data::Dumper;
use File::Copy;

use vars qw/$opt_b $opt_c $opt_m $opt_a $opt_p/;

getopts("b:c:m:a:p:");

my $usage = "Usage:  $0 \
                -b list of blast output files\
                -c list of clusters from uclust (optional)\ 
                -m input mapping file\
		-a annotation file\
                -p output prefix for the processed files\
                \n";

die $usage unless defined $opt_b
              and defined $opt_m
              and defined $opt_a
              and defined $opt_p;
#
my $list           = $opt_b;
my $mapfile        = $opt_m;
my $annotationfile = $opt_a;
my $prefix         = $opt_p;
my $clusterlist    = "";
if (defined($opt_c)){
  $clusterlist = $opt_c;
}
#
my %samples         = ();
my %mapdata         = ();
my %map_levels      = ();
my %map_types       = (); 
my $num_map_levels  = 0;
my %antndata        = ();
my %antn_levels     = ();
my %antn_types      = ();
my $num_antn_levels = 0;
#
my @orderedsamples  = (); 
my %clusters        = ();
my %COUNTS          = ();
#

# load up the mapping_data
open IN, "$mapfile" or die "Can't open mapping file: $mapfile\n";
while(<IN>){
  chomp($_);
  next if ($_ eq "");
  my @A = split "\t", $_;
  if ($_ =~ /^\#File/){
    $num_map_levels = $#A+1-1;            
    for my $i (1 .. $#A){
      $map_levels{$i} = $A[$i];
    }
  }else{
    $samples{$A[0]} = 0;
    for my $i (1 .. $#A){
      $mapdata{$A[0]}{$i} = $A[$i];
      $map_types{$i}{$A[$i]} = 0;   
    }
  }  
}
close IN;


# order samples by the first mapping feature
if ($num_map_levels == 0){
  @orderedsamples = sort keys %samples;  
}else{
  foreach my $t (sort keys %{$map_types{1}}){
    foreach my $s (keys %samples){
      if ($mapdata{$s}{1} eq $t){
        push @orderedsamples, $s; 
      }
    }
  }    
}


# load up annotation data
open IN, "$annotationfile" or die "Can't open annotation file: $annotationfile\n";
while(<IN>){
  chomp($_);
  next if ($_ eq "");
  my @A = split "\t", $_;
  if ($_ =~ /^Name/){
    $num_antn_levels = $#A+1-1;
    for my $i (1 .. $#A){
      $antn_levels{$i} = $A[$i];
    }
  }else{
    for my $i (1 .. $#A){
      $antndata{$A[0]}{$i} = $A[$i];
    }
  }
}
close IN;


# load up cluster info if necessary
if ($clusterlist ne ""){
  my $catclustlist = `cat $clusterlist`;
  my @catclustlist = split "\n", $catclustlist;
  for my $i (0 .. $#catclustlist){
    open IN, "$catclustlist[$i]" or die "Can't open $catclustlist[$i] for processing!\n"; 
    my $seed = "";
    while(<IN>){
      chomp($_);
      next if ($_ eq "");
      next if ($_ =~ /^>Cluster/);
      my @A = split ">", $_;
      my @B = split /\.\.\./, $A[1]; 
      if ($_ =~ /^0/){
        $seed = $B[0];
      }else{
        push @{$clusters{$seed}}, $B[0];
      }
    }
  }
}

# how many blast files are provided?
my $listlength = `wc $list`;
my @listlength = split " ", $listlength;

# one hit per seq
my %query_hit = ();

# we assume that each file represents a different specific sample
my $catstr = `cat $list`;
my @catstr = split "\n", $catstr;
for my $i (0 .. ($listlength[0]-1)){
  # do some processing of the filename to get the prefix
  # and store the associated barcode
  my @line = split /\//, $catstr[$i];
  my $filename = $line[$#line];

  #open and process this file
  open IN, "$catstr[$i]" or die "Can't open $catstr[$i] for preprocessing!\n";
  while(<IN>){
    chomp($_);
    my @blastline = split "\t", $_;

    next if (defined($query_hit{$blastline[0]})); # one hit allowed per sequence
    $query_hit{$blastline[0]} = 1; # if you havent seen this query before, catalog it

    my @qname = split /\_/, $blastline[0];
    my $qname = join("", @qname[0..($#qname-1)]);     
    my $hit   = $blastline[1];
    
    if (!defined($antndata{$hit})){ # if we don't have a record for it, don't count it
      next;
    }

    # catalog the hit in the COUNTS datatype
    foreach my $a (keys %antn_levels){
      $COUNTS{$qname}{$a}{$antndata{$hit}{$a}}++;   
      $antn_types{$a}{$antndata{$hit}{$a}} = 1; 

      # also count this hit for any seqs in the associated cluster
      if (defined($clusters{$blastline[0]})){
        foreach my $qclust (@{$clusters{$blastline[0]}}){
          my @sampstr = split /\_/, $qclust;
          my $sampstr  = join("", @sampstr[0..($#sampstr-1)]);
          if (!defined($COUNTS{$sampstr}{$a}{$antndata{$hit}{$a}})){
            $COUNTS{$sampstr}{$a}{$antndata{$hit}{$a}} = 1;
          }else{
            $COUNTS{$sampstr}{$a}{$antndata{$hit}{$a}}++;  
          }
        }
      }
    }
  }
  close IN;
  # end of this file
  # move on to the next file
} 
# End of all blast list


# zero out empty entries
foreach my $a (keys %antn_levels){
foreach my $s (keys %samples){
foreach my $t (keys %{$antn_types{$a}}){
  if (!defined($COUNTS{$s}{$a}{$t})){
    $COUNTS{$s}{$a}{$t} = 0;
  }
}  
}
}

#
# now print out the revelant tables
foreach my $a (keys %antn_levels){
  printTable($a);
}

# now if there are map_types that require pairwise comparisons
# print out those tables now
foreach my $aa (keys %antn_levels){
foreach my $a (keys %map_levels){
  next if (substr($map_levels{$a},-2) ne "_p"); #if this -- then they don't want a pairwise comparison for this

  my @gs = sort keys %{$map_types{$a}};
  for my $i (0 .. $#gs){
    for my $j (0 .. $#gs){
      next if ($j >= $i);
      printPairedGroups($a, $aa, $gs[$i], $gs[$j]);
    }
  }
}  
}

#**********************************************************************************
sub printPairedGroups
{
  my ($maplevel, $antnlevel, $g1, $g2) = @_;

  my @pairedorderedsamples = ();
  my $g1count = 0;
  for my $s (0 .. $#orderedsamples){
    if ($mapdata{$orderedsamples[$s]}{$maplevel} eq $g1){
      push @pairedorderedsamples, $orderedsamples[$s];
      $g1count++; 
    }    
  }

  my $g2count = 0;
  for my $s (0 .. $#orderedsamples){
    if ($mapdata{$orderedsamples[$s]}{$maplevel} eq $g2){
      push @pairedorderedsamples, $orderedsamples[$s];
      $g2count++;
    }
  }

  return if ($g1count <= 0 or $g2count <=0);
  return if (($g1count == 1 and $g2count != 1) or ($g1count != 1 and $g2count == 1));  


  open OUT, ">$prefix.$antn_levels{$antnlevel}.$g1\_vs_$g2.$g1count-$g2count.2tsv" or die "Can't open $prefix.$antn_levels{$antnlevel}.$g1\_vs_$g2.$g1count-$g2count.2tsv for writing!\n";
  for my $s (0 .. $#pairedorderedsamples){
    print OUT "\t$pairedorderedsamples[$s]";
  }
  print OUT "\n";

  # print feature tables 
  foreach my $f (sort keys %{$antn_types{$antnlevel}}){
    # if it's all zeros, then skip it:
    my $fsum = 0;
    for my $i (0 .. $#pairedorderedsamples){
      $fsum += $COUNTS{$pairedorderedsamples[$i]}{$antnlevel}{$f};
    } 
    next if ($fsum == 0);
  
    # if we get to here, we can print it out ~
    print OUT "$f";
    for my $i (0 .. $#pairedorderedsamples){
      print OUT "\t$COUNTS{$pairedorderedsamples[$i]}{$antnlevel}{$f}";
    }
    print OUT "\n";
  }
  close OUT;
}


sub printTable
{
  my ($level) = @_;
  # open the file
  open OUT, ">$prefix.$antn_levels{$level}.tsv" or die "Cannot open $prefix.$antn_levels{$level}.tsv for writing...\n";

  # print header line
  for my $i (0 .. $#orderedsamples){
    print OUT "\t$orderedsamples[$i]";
  }  
  print OUT "\n";
 
  # print feature tables 
  foreach my $f (sort keys %{$antn_types{$level}}){
    print OUT "$f";
  for my $i (0 .. $#orderedsamples){
    print OUT "\t$COUNTS{$orderedsamples[$i]}{$level}{$f}";
  }
    print OUT "\n";
  }
  close OUT;
}





