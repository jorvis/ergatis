#!/usr/bin/env perl
use strict;
use warnings;

#******************************************************************************
# *clovr_metagenomics_tsvtables.pl
# Author: james robert white, james.dna.white@gmail.com

# This function takes a list blast outputs (-m 8) and some metadata and creates
# a set of files for postprocessing in the clovr metagenomics pipeline.

# A list of 1 or more blast output files is given
# Samples are described by the associated mapping file. Hits are described by
# the annotation file provided as input.
# 
# The outputs of this program are tab-delimited tables of hits
# Also outputs internal read names + annotations in tab delimited form
# for cross-annotation later on
#******************************************************************************
use Getopt::Std;
use Data::Dumper;
use File::Copy;

use vars qw/$opt_b $opt_c $opt_m $opt_a $opt_f $opt_y $opt_z $opt_p/;

getopts("b:c:m:a:f:y:z:p:");

my $usage = "Usage:  $0 \
                -b list of blast output files\
                -c list of clusters from uclust (optional for clovr-metagenomics)\ 
                -m input mapping file\
		-a annotation file\
                -f fasta file of blasted sequences (optional - for clovr-its)\
                -y bsml list from metagene runs (optional)\
                -z list of polypeptide clusters from uclust (optional)
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
  $clusterlist     = $opt_c;
}

my $bsml_list      = "";
if (defined($opt_y) and $opt_y ne ""){
  $bsml_list       = $opt_y;
}

my $polypep_clusterlist = "";
if (defined($opt_z) and $opt_z ne ""){
  $polypep_clusterlist  = $opt_z;
}

my $fasta_file = ""; # of blasted sequences
if (defined($opt_f) and $opt_f ne ""){
  $fasta_file  = $opt_f;
}

#
my %samples         = ();
my %mapdata         = ();
my %map_levels      = ();
my %map_types       = (); 
my $num_map_levels  = 0;
my %antndata        = ();
my %antn_levels     = ();
my @ordered_antn_levels = ();
my %antn_types      = ();
my $num_antn_levels = 0;
my %polyrepmap      = ();  # maps translated polypeptides to orig reads
my %polyclusters    = ();
#
my @orderedsamples  = (); 
my %clusters        = ();
my %COUNTS          = ();
my $QIIMEOTUs       = 0;  # flag for additional processing for CloVR_ITS
my %lengths         = (); # if fasta file is provided, this will store the lengths of the sequences
#
if ($bsml_list eq "" and $polypep_clusterlist ne ""){
  die "**ERROR** a polypeptide clustering is provided but no bsml list of initial polypeptides!!\n";
}
#

#*****************************************************************
# ************* DATA LOADING ************************************#
#*****************************************************************

# load of polypeptide mapping to DNA sequences
if ($bsml_list ne ""){
  print "loading polypeptide map\n";
  loadUpPolyMap($bsml_list);
}

# load up polypeptide clusters from UCLUST
if ($polypep_clusterlist ne ""){
  print "loading polypeptide clusters\n";
  loadPolyPepClusters($polypep_clusterlist);    
}

# load sequence lengths
if ($fasta_file ne ""){
  print "loading up sequence lengths\n";
  open IN, "$fasta_file" or die "Can't open $fasta_file\n";
  my $cseq  = "";
  my $cname = "";
  while(<IN>){
    chomp($_);
    if ($_ =~ />/){
      if ($cseq ne ""){
        $lengths{$cname} = length($cseq);
      }
      my @A = split " ", $_;
      $cname = substr($A[0], 1);
      $cseq  = "";
    }else{
      $cseq .= $_;   
    }
  }
  close IN;
  # last sequence in the file!
  $lengths{$cname} = length($cseq);
}

# load up cluster information
if ($clusterlist ne "" and $clusterlist ne "none"){
  loadUpDNAClusters();
}

# load up the mapping_data
loadMappingData();

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

# load up annotation data from the db
loadAnnotationData();

# define the raw output annotation to be written
open RAW, ">$prefix.raw.antn" or die "Can't open $prefix.raw.antn for writing!\n";
print RAW "SequenceName";
foreach my $aa (@ordered_antn_levels){
  print RAW "\t$aa";
}
print RAW "\n";

#*******************************************************************
#***** BEGIN BLAST FILE ANALYSIS **********************************#
#*******************************************************************
# how many blast files are provided?
my $listlength = `wc $list`;
my @listlength = split " ", $listlength;

# keep track -- one hit per blasted seq
# if the sequence is already defined in this
# hash it's already been recorded...
my %query_hit = ();

# recall that $list is the list of raw blast output files
my $catstr = `cat $list`;
my @catstr = split "\n", $catstr;
for my $i (0 .. ($listlength[0]-1)){

  # get the filename within the directory
  my @line = split /\//, $catstr[$i];
  my $filename = $line[$#line];

  # open and process this file
  open IN, "$catstr[$i]" or die "Can't open $catstr[$i] for preprocessing!\n";
  while(<IN>){
    chomp($_);
    my @blastline = split "\t", $_;

    next if (defined($query_hit{$blastline[0]})); # one hit allowed per sequence
    $query_hit{$blastline[0]} = 1;                # if you havent seen this query before, catalog it, we output multiple hits now

    # scenario 1:
    if ($bsml_list eq "" and $polypep_clusterlist eq ""){ # then there are no clusters to think about
      catalog($blastline[0], $blastline[1], $blastline[2], $blastline[3]);
    }

    # scenario 2:
    # if just a bsml file list is provided then convert the hit name to a rep DNA sequence:
    if ($bsml_list ne "" and $polypep_clusterlist eq ""){
      catalog($polyrepmap{$blastline[0]}, $blastline[1], $blastline[2], $blastline[3]); 
    }
  
    # scenario 3:
    # if a bsml file list is provided and a polypep cluster list then account
    # for all propagated hits: a peptide spreads to the whole cluster and each
    # peptide in the cluster goes to a rep sequence...
    if ($bsml_list ne "" and $polypep_clusterlist ne ""){
      catalog($polyrepmap{$blastline[0]}, $blastline[1], $blastline[2], $blastline[3]);  
      foreach my $pclustmember (@{$polyclusters{$blastline[0]}}){
        catalog($polyrepmap{$pclustmember}, $blastline[1], $blastline[2], $blastline[3]);  
      }
    }
  }
  close IN;

  # end of this file
  # move on to the next file
}

# finally if this is CloVR-ITS, we also want to report any sequences that had no BLAST hits at all as unclassified
if ($QIIMEOTUs == 1){
  foreach my $s (sort keys %lengths){
    if (!defined($query_hit{$s})){
      print "$s had no blast hit...\n";
      catalog($s, 212374277, 0, 0);
    }
  }
}
 
close RAW;

#***************************************************************
#***************************** End of all blast files *********#
#***************************************************************

zeroOutEmptyEntries();

# now print out the revelant tables
foreach my $a (keys %antn_levels){
  printTable($a);
}

#***********************************************************************************
# PAIRED GROUP TABLES
#***********************************************************************************
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
# Subroutines
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

# takes a list of BSML files and 
# creates a mapping from DNA reads to predicted polypeptides from Metagene
sub loadUpPolyMap
{
  my ($bsmlliststr) = @_;
  my $catbsmllist = `cat $bsmlliststr`;
  my @catbsmllist = split "\n", $catbsmllist;
  for my $i (0 .. $#catbsmllist){
    open IN, "$catbsmllist[$i]" or die "Can't open $catbsmllist[$i] for processing!\n";
    my $rep = "";
    while(<IN>){
      chomp($_);
      next if ($_ eq "");

      if ($_ =~ /<Sequence length=/){
        #  <Sequence length="999" class="assembly" id="sludgeUS.small.fna_49" molecule="dna">
        my @A = split /id=\"/, $_;
        my @B = split /\"/, $A[1];
        $rep = $B[0];
      }

      if ($_ =~ /<Feature class="polypeptide" id=/){
        #  <Feature class="polypeptide" id="new.polypeptide.X88029X879.1">
        my @A = split /id=\"/, $_;
        my @B = split /\"/, $A[1];
        my $poly = $B[0];
        $polyrepmap{$poly} = $rep;         
      }
    }
  }
}

# takes a list of cluster files
# from UCLUST and defined the polyclusters type
# a peptide points to an array of others representing
# the cluster
sub loadPolyPepClusters
{
  my ($pp_clusterlist) = @_;

  my $catclustlist = `cat $pp_clusterlist`;
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
        push @{$polyclusters{$seed}}, $B[0];
      }
    }
  }
}


sub catalog
{
  my ($repseqhit, $hit, $alignid, $alignlength) = @_;

  my @qname = split /\_/, $repseqhit;
  my $qname = join("_", @qname[0..($#qname-1)]);

  my %clovr_its_thresholds = ("Class"   => 60,   
                              "Order"   => 70,   
                              "Family"  => 75,  
                              "Genus"   => 85,   
                              "Species" => 90);

  if (!defined($antndata{$hit})){ # if we don't have a record for it, then the annotation doesn't cover it, so don't count it
    next;
  }

  my %mini_antn_hash = (); # This will hold the annotations of each sequence to be printed out later
  my $a = 0;

  # catalog the hit in the COUNTS datatype
  # for each level in the annotation hierarchy (e.g. Phylum, Class, etc)
  foreach my $levelname (@ordered_antn_levels){ 
    $a++;
    if (defined($mini_antn_hash{$repseqhit})){ # then it's time to add a tab to the string
      $mini_antn_hash{$repseqhit} .= "\t";     # tab-delimit between levels
  
      if (defined($clusters{$repseqhit})){ # then the representative seq represents a cluster of sequences
        foreach my $qclust (@{$clusters{$repseqhit}}){
          $mini_antn_hash{$qclust} .= "\t";   
        }
      }
    } 

  # this ref sequence may have multiple annotations at this level, so all will be added
  foreach my $b (@{$antndata{$hit}{$a}}){ 
    my $annotationToAdd = $b;
 
    if ($QIIMEOTUs == 1){ # then this is a clovr-its run, and we need to decide if it meets the threshold criteria
      my $repseqcoverage = 100*$alignlength/$lengths{$repseqhit}; 
      if ($repseqcoverage < 90 or $alignid < $clovr_its_thresholds{$levelname}){
        $annotationToAdd = "Unclassified";
      } 
    }

    $COUNTS{$qname}{$a}{$annotationToAdd}++;
    $antn_types{$a}{$annotationToAdd} = 1;

    if (!defined($mini_antn_hash{$repseqhit})){ # then start the string
      $mini_antn_hash{$repseqhit}  = $annotationToAdd;
    }else{ # the string has been started
      my $lastchar = substr($mini_antn_hash{$repseqhit},-1);
      if ($lastchar eq "\t"){ # then this is a new level of annotation
        $mini_antn_hash{$repseqhit} .= $annotationToAdd;
      }else{ # then it's not a new level of annotation
        $mini_antn_hash{$repseqhit} .= ",$annotationToAdd";
      }
    }

    # this sequence hit, $repseqhit may be associated with a cluster of DNA seqs
    # so, also count this hit for any seqs in the associated nucleotide sequence cluster
    if (defined($clusters{$repseqhit})){
      foreach my $qclust (@{$clusters{$repseqhit}}){

        if (!defined($mini_antn_hash{$qclust})){
          $mini_antn_hash{$qclust}  = $annotationToAdd;
        }else{
          my $lastchar = substr($mini_antn_hash{$qclust},-1);
          if ($lastchar eq "\t"){ # then this is a new level of annotation
            $mini_antn_hash{$qclust} .= $annotationToAdd;
          }else{ # then it's not a new level of annotation
            $mini_antn_hash{$qclust} .= ",$annotationToAdd";
          }
        }

        my @sampstr  = split /\_/, $qclust;
        my $sampstr  = join("", @sampstr[0..($#sampstr-1)]);
        if (!defined($COUNTS{$sampstr}{$a}{$annotationToAdd})){
          $COUNTS{$sampstr}{$a}{$annotationToAdd} = 1;
        }else{
          $COUNTS{$sampstr}{$a}{$annotationToAdd}++;
        }
      }
    }
  }
  }

  # print out the mini hash 
  foreach my $mem (keys %mini_antn_hash){
    print RAW "$mem\t$mini_antn_hash{$mem}\n";
  }

}


sub loadMappingData
{
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
      if ($QIIMEOTUs == 1 or $opt_c eq "none" or !defined($opt_c)){ # then this is for clovr_its with no clusters, so just keep the prefixes
        my @B = split /\./, $A[0];
        $A[0] = join(".", @B[0 .. ($#B-1)]);  
      }
      $samples{$A[0]} = 0;
      for my $i (1 .. $#A){
        $mapdata{$A[0]}{$i} = $A[$i];
        $map_types{$i}{$A[$i]} = 0;
      }
    }
  }
  close IN;

}


sub loadAnnotationData
{
  open IN, "$annotationfile" or die "Can't open annotation file: $annotationfile\n";
  while(<IN>){
    chomp($_);
    next if ($_ eq "");
    my @A = split "\t", $_;
    if ($_ =~ /^Name/){
      $num_antn_levels = $#A+1-1;
      for my $i (1 .. $#A){
        $antn_levels{$i} = $A[$i];
        push @ordered_antn_levels, $A[$i];
      }
    }else{
      for my $i (1 .. $#A){
        my @B = split /\$/, $A[$i]; # if there is more than one annotation at this level, 
                                    #  then we delimit by $ symbols.
        for my $b (0 .. $#B){
          push @{$antndata{$A[0]}{$i}}, $B[$b];
        }
      }
    }
  }
  close IN;
}



sub zeroOutEmptyEntries
{
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
}

# defines cluster type
# a rep dna sequence points to a cluster of other dna sequences
sub loadUpDNAClusters
{
  my $catclustlist = `cat $clusterlist`;
  my @catclustlist = split "\n", $catclustlist;
  for my $i (0 .. $#catclustlist){

    my $toplineck = `head -n 1 $catclustlist[$i]`;
    if ($toplineck =~ /^0\t/){
      $QIIMEOTUs = 1;
    } 

    open IN, "$catclustlist[$i]" or die "Can't open $catclustlist[$i] for processing!\n";
    my $seed = "";
    while(<IN>){
      chomp($_);
      next if ($_ eq "");

      if ($QIIMEOTUs == 0){

        next if ($_ =~ /^>Cluster/);
        my @A = split ">", $_;
        my @B = split /\.\.\./, $A[1];
        if ($_ =~ /^0/){
          $seed = $B[0];
        }else{
          push @{$clusters{$seed}}, $B[0];
        }
      }else{ # It's a QIIME based OTU file from CloVRITS
        my @A = split "\t", $_;
        for my $j (1 .. $#A){
          push @{$clusters{$A[0]}}, $A[$j];
        }
      }  
    }
  }

}

