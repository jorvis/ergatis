#! /usr/local/bin/perl -w


use strict;

use lib "/home/jravel/lib/";
use BeginPerlBioinfo;

my $dir4 = $ARGV[0];

unless(opendir(DIRECTORY, "$dir4")) {   
    print "Cannot open directory $dir4";
    exit;
}

my @dirs = grep (!/^\.\.?$/, readdir(DIRECTORY));  

closedir(DIRECTORY);

foreach my $dirlist (@dirs) {

    chomp $dirlist;

    if ($dirlist =~ /^GBA/) {

	my ($tag1, $pre) = ($dirlist =~ /^(GBA\d*)_\D{3}_\d_\d*(\.\w*)/); 

	my $newname = $tag1 . $pre;
	system ("mv $dir4/$dirlist $dir4/$newname");
    }else {
	next;
    }
}

exit;

