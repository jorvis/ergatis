#!/usr/bin/perl -w

use strict;

my $file = shift || die "pass me a file";

open my $ifh, $file || die "can't open file $file : $!";

my $in_dcespec = 0;

while (<$ifh>) {
    
    ## do we match an open tag?
    if ( m|\<dceSpec| ) {
        $in_dcespec = 1;
    }
    
    if (! $in_dcespec) {
        print;
    }

    ## do we match a close tag?
    if ( m|\</dceSpec| ) {
        $in_dcespec = 0;
    }        
}
