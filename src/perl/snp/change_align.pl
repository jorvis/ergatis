#!/usr/local/bin/perl -w

use strict;
use lib "/home/jravel/lib/";
use BeginPerlBioinfo;

my @F;
my $infile = $ARGV[0];

my @file = get_file_data ( $infile);
my $marker = 0;
my $i =0;
foreach my $line (@file) {
    
#    print $line;

    if ($line =~ /^>/ && $marker == 0) {

	$marker = 1;
	next;
    
    }elsif ( $line =~ /Ref/ && $marker == 1) {

	@F = split(" ", $line);

	print "$F[2]\t$F[4]\t$F[5]\t$F[3]\n";

	$i++;

	next;

    }elsif ($line =~ /^>/ && $marker == 1) {

	last;
    }
}


print "Total number of errors = $i\n";

exit;



    
