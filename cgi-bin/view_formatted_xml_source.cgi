#!/local/perl/bin/perl -w

use strict;
use CGI;
use CGI::Carp qw(fatalsToBrowser);

my $q = new CGI;

print $q->header( -type => 'text/html' );

## will be like:
## /usr/local/scratch/annotation/TGA1/Workflow/split_fasta/29134_test2/pipeline.xml
my $file = $q->param("file") || die "pass file";

pageHeader();

## don't do it if the file doesn't end in .xml or .instance
if ($file !~ /\.xml$/ && $file !~ /\.instance$/ && $file !~ /\.bsml$/) {
    print STDERR "skipped display of $file in source viewer\n";
    quitNicely("i decline to show this type of file.");
}

print "<pre>";

## open the file and print it to the screen.
open (my $ifh, "<$file") || quitNicely("couldn't open file $file");

while (my $line = readline $ifh) {
    my ($tag, $url, $word);
    
    ## colorize any tags
    my $newline = $line;
    while ( $line =~ m|\<(/*\w+)|g ) {
        $tag = $1;
        $newline =~ s|$tag|STARTTAGSPAN${tag}ENDTAGSPAN|;
    }

    $line = $newline;

    ## look for comments
    $line =~ s/<!--/STARTCOMMENTSPAN<!--/g;
    $line =~ s|-->|-->ENDCOMMENTSPAN|g;
    
    $line =~ s/\</\&lt\;/g;
    $line =~ s/\>/\&gt\;/g;

    ## add the spans properly
    $line =~ s|STARTTAGSPAN|<span class="tag">|g;
    $line =~ s|STARTCOMMENTSPAN|<span class="comment">|g;
    $line =~ s!(ENDTAGSPAN|ENDCOMMENTSPAN)!</span>!g;

    
    ## look for any linkable xml
    if ( $line =~ m^(?<!\$\;)(/[/a-z0-9_\-.]+\.(?:xml|instance|bsml))(?![\./])^i ) {
        $url = $1;
        $line =~ s|$url|<a href="./view_formatted_xml_source.cgi?file=$url">$url</a>|;
    }
    
    ## look for any linkable ini
    if ( $line =~ m^(?<!\$\;)(/[/a-z0-9_\-.]+\.(?:ini|config|conf))(?![\./])^i ) {
        $url = $1;
        $line =~ s|$url|<a href="./view_formatted_ini_source.cgi?file=$url">$url</a>|;
    }

    ## look for any linkable lists
    if ( $line =~ m^(?<!\$\;)(/[/a-z0-9_\-.]+\.list)(?![\./])^i ) {
        $url = $1;
        $line =~ s|$url|<a href="./view_raw_source.cgi?file=$url">$url</a>|;
    }
    
    print $line;
}

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
        span.tag {
            color: rgb(100,100,100);
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
