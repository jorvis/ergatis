#! /usr/local/bin/perl -w


use strict;
use lib "/home/jravel/lib";
use BeginPerlBioinfo;


my $pos = $ARGV[0];

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

print  $pos, "\n";

exit;
