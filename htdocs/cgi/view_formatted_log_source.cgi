#!/usr/bin/perl -w

use strict;
use CGI;
use CGI::Carp qw(fatalsToBrowser);

my $q = new CGI;

print $q->header( -type => 'text/html' );

## will be like:
## /usr/local/scratch/annotation/TTA1/workflow/runtime/wu-blastp/20724_AllGroup.niaa/component.conf.bld.ini
my $file = $q->param("file") || die "pass file";

pageHeader();

## don't do it if the file doesn't end in .ini or .conf or config
if ($file !~ /\.out$/ && $file !~ /\.log$/ && $file !~ /\.stderr$/ && $file !~ /\.stdout$/) {
    print STDERR "skipped display of $file in source viewer\n";
    quitNicely("i decline to show this type of file.");
}

print "<pre>";

## open the file and print it to the screen.
open (my $ifh, "<$file") || quitNicely("couldn't open file $file");

while (my $line = readline $ifh) {
    chomp $line;
    my ($url, $var);

    ## colorize any comments
    if ($line =~ /^[\#\;]/) {
        $line = '<span class="comment">' . $line . '</span>';
    }
    
    ## colorize variables
    my $newline = $line;
    while ($line =~ /\$\;(\w+)\$\;/g) {
        $var = $1;
        $newline =~ s|\$\;$var\$\;|\$\;<span class="inivar">$var</span>\$\;|;
    }
    $line = $newline;
    
    
    ## embolden section headers
    if ($line =~ /^\[/) {
        $line = '<span class="inihead">' . $line . '</span>';
    }

    ## look for any linkable xml
    if ( $line =~ m^(?<!\$\;)(/[/a-z0-9_\-.]+\.(?:xml|instance|bsml))\s*$^i ) {
        $url = $1;
        $line =~ s|$url|<a href="./view_formatted_xml_source.cgi?file=$url">$url</a>|;
    }
    
    ## look for any linkable ini
    if ( $line =~ m^(?<!\$\;)(/[/a-z0-9_\-.]+\.(?:ini|config|conf))\s*$^i ) {
        $url = $1;
        $line =~ s|$url|<a href="./view_formatted_ini_source.cgi?file=$url">$url</a>|;
    }
    
    ## look for any linkable lists
    if ( $line =~ m^(?<!\$\;)(/[/a-z0-9_\-.]+\.list)\s*$^i ) {
        $url = $1;
        $line =~ s|$url|<a href="./view_raw_source.cgi?file=$url">$url</a>|;
    }
    
    print "$line\n";
}

my %h;
map { $h{$_} = 1 } ('perl.c', 'sv.c', 'hv.c', 'av.c');
for ('perl.c', 'sv.c', 'hv.c', 'av.c') { $h{$_} = 1 } 


print "</pre>";
print "</body></html>";

sub pageHeader {
    print <<heADerStuff;

<html>
<head>
    <style type="text/css">
        #file {
            background-color: rgb(0,0,150);
            color: rgb(255,255,255);
            font-family: verdana, helvetica, arial, sans-serif;
            font-size: 10px;
            font-weight: bold;
            padding: 5px 0px 5px 10px;
        }
        span.comment {
            background-color: rgb(200,200,200);
        }
        span.inihead {
            font-weight: bold;
        }
        span.inivar {
            color: rgb(0,75,0);
        }
        a {
            font-weight: bold;
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
