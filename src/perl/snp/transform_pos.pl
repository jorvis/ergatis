#! /usr/local/bin/perl -w

### THIS SCRIPT WILL TRANSFORM A COORDINATE FROM 
### THE 6615 assembly into a 6609 assembly of B. anthracis Ames
###
### 6615 is the annotation molecule 
### 6609 is the closure molecule

use strict;
use lib "/home/jravel/lib";
use BeginPerlBioinfo;


my $infile = $ARGV[0];

my @data = get_file_data($infile);

open (OUT, ">mum_trans.out");


foreach my $line (@data) {

    my @temp = split("\t", $line);




    my $pos = $temp[0];

    if ($pos =~ /^Total/) {
	last;
    }
    
    if ($pos < 1013547) {
	$pos = $pos + 576235;
    }elsif ( $pos > 1013546 && $pos < 1288255) {
	$pos = $pos + 576234;
    }elsif ( $pos > 1288254 && $pos < 3390649) {
	$pos = $pos + 576235;
    }elsif ( $pos > 3390648 && $pos < 4651059) {
	$pos = $pos + 576236;
    }elsif ( $pos > 4651058) {
	$pos = $pos - 4651058;
    }




    my $newline = "$pos\t$temp[1]\t$temp[2]\t$temp[3]";
    
    print OUT $newline;
    
}


close (OUT);

exit;
