#!/usr/bin/perl -w

use strict;
use CGI;
use CGI::Carp qw(fatalsToBrowser);
use Ergatis::ConfigFile;

my $q = new CGI;

print $q->header( -type => 'text/plain' );

## will be like:
my $file = $q->param("file") || die "pass file";


## the file may have been compressed
if ( ! -e $file && -e "$file.gz" ) {
    $file .= '.gz';
}

## read the ergatis config file
my $ergatis_cfg = new Ergatis::ConfigFile( -file => "ergatis.ini" );

## only print certain file types (for minor security)
if ($file !~ /\.xml$/ && 
    $file !~ /\.instance$/ && 
    $file !~ /\.bsml$/ &&
    $file !~ /\.gz$/ && 
    $file !~ /\.ini/ &&
    $file !~ /\.log/ &&
    $file !~ /\.config/ &&
    $file !~ /\.list/ &&
    $file !~ /\.stderr/ &&
    $file !~ /\.stdout/ && 
    ! $ergatis_cfg->val('display_settings', 'show_all_files')) {
    print STDERR "skipped display of $file in source viewer\n";
    quitNicely("i decline to show this type of file.");
}

## open the file and print it to the screen.
my $ifh;
if ($file =~ /\.gz$/) {
    open($ifh, "<:gzip", $file) || quitNicely("couldn't open compressed file $file: $!");
} else {
    open($ifh, "<$file") || quitNicely("couldn't open file $file");
}

while (my $line = readline $ifh) {
    ## just print it    
    print $line;
}

sub quitNicely {
    my $msg = shift;
    
    print $msg;
    exit;
}
