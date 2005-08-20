#!/usr/local/bin/perl -w

use strict;
use CGI;
use CGI::Carp qw(fatalsToBrowser);

my $q = new CGI;

print $q->header( -type => 'text/html' );

## will be like:
my $file = $q->param("file") || die "pass file";

pageHeader();

## only print certain file types (for minor security)
if ($file !~ /\.list$/) {
    print STDERR "skipped display of $file in source viewer\n";
    quitNicely("i decline to show this type of file.");
}

print "<pre>";

## open the file and print it to the screen.
open (my $ifh, "<$file") || quitNicely("couldn't open file $file");

while (my $line = readline $ifh) {
    ## just print it    
    print $line;
}

print "</pre>";
print "</body></html>";

sub pageHeader {
    print <<heADerStuff;

<html>
<head>
    <style type="text/css">
        body {
            font-family: verdana, helvetica, arial, sans-serif;
            font-size: 10px;
        }
        #file {
            background-color: rgb(0,0,150);
            color: rgb(255,255,255);
            font-weight: bold;
            padding: 5px 0px 5px 10px;
        }
    </style>
</head>
<body>

<div id="file">
    $file
</div>

heADerStuff
}

sub quitNicely {
    my $msg = shift;
    
    print $msg;
    exit;
}
