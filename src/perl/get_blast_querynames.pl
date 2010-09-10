#!/usr/bin/perl

eval 'exec /usr/bin/perl  -S $0 ${1+"$@"}'
    if 0; # not running under some shell
use strict;
use warnings;

#******************************************************************************
# *get_blast_querynames.pl
# Author: james robert white, james.dna.white@gmail.com

# This function takes a list of blast files as input.
# The blast files msut be in -m 8 output. This script only grabs the names.
#******************************************************************************
use Getopt::Std;
use Data::Dumper;
use warnings;

use vars qw/$opt_f $opt_o/;

getopts("f:o:");

my $usage = "Usage:  $0 \
                -f <list of blast output files>\
                -o <output file name>\
                \n";

die $usage unless defined $opt_f
              and defined $opt_o;

my $outputfile = $opt_o;
my $blastlist  = $opt_f;
my @blastfiles = ();

# now open up the fasta list and process the files
if (defined($blastlist)) {
  open (FLIST, $blastlist) or die ("Could not open FASTA list $blastlist: $!");
  while (my $file = <FLIST>) {
    chomp($file);
    push (@blastfiles, $file);
  }
  close (FLIST);
}

foreach my $file (@blastfiles){  
  `cut -f 1 $file >> $outputfile`; 
}





