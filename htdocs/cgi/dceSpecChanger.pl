#!/usr/local/bin/perl -w

use strict;

my $file = shift || die "pass me a file";

my $ifh;
my $ofh;
if( $file =~ /\.gz$/ ) {
    open $ifh, "<:gzip", $file || die "can't open file $file : $!";
    open $ofh, ">:gzip", "$file.tmp" || die "can't open $file.tmp for writing : $!";
} else {
    open $ifh, $file || die "Can't open file $file : $!";
    open $ofh, "> $file.tmp" || die "can't open $file.tmp for writing : $!";        
}

my ($project_code, $passthrough);

my $in_dcespec = 0;

while (<$ifh>) {
    
    ## do we match an open tag?
    if ( m|\<dceSpec| ) {
        $in_dcespec = 1;
    }
    
    if (! $in_dcespec) {
        print $ofh $_;
    } elsif ( /\<group\>([^\<]+)\</ ) {
        $project_code = $1;
    } elsif ( /\<passthrough\>([^\<]+)\</ ) {
	$passthrough = $1;
    }

    ## do we match a close tag?
    if ( m|\</dceSpec| ) {
        $in_dcespec = 0;
        print $ofh "<dceSpec type=\"sge\"><OS>linux</OS><group>$project_code</group><passthrough>$passthrough</passthrough></dceSpec>\n";
    }        
}

system("mv $file.tmp $file");
close($ofh);
