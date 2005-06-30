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

## give colors as rgb values or hexidecimal
my %colors = (
                complete    => 'rgb(0,200,0)',      ## green
                incomplete  => 'rgb(75,75,75)',     ## dark grey
                failed      => 'rgb(200,0,0)',      ## red
                pending     => 'rgb(200,200,0)',    ## yellow
                errors      => 'rgb(200,0,0)',      ## red
                error       => 'rgb(200,0,0)',      ## red
                running     => 'rgb(0,0,200)',      ## blue
                waiting     => 'rgb(200,200,0)',    ## yellow
                interrupted => 'rgb(200,0,200)',    ## purple
                total       => 'rgb(0,0,0)',        ## black 
             );
my $progress_image_width = 500;

print "<div id='sourcecode'>";

## open the file and print it to the screen.
open (my $ifh, "<$file") || quitNicely("couldn't open file $file");

## something to remember what states we find
my %states;
my $overall_state = 0;
my $found_states = 0;
my $within_status_box = 0;
my $command_count = 0;

while (my $line = readline $ifh) {
    my ($tag, $url, $word);

    ## have we found the first state?
    if ( (! $overall_state) && $line =~ m|<state>(.+)</state>|) {
        $overall_state = $1;
    }

    ## are we dealing with status box?
    if ($line =~ /<status>/) {
        $within_status_box = 1;
    } elsif ($line =~ /<\/status>/) {
        $within_status_box = 0;
        $found_states = 1;
    } elsif ($within_status_box && (! $found_states) && $line =~ /\<(.+?)\>(\d+)/) {
        if ($2) {
            $states{$1} += $2;
            $command_count += $2 unless ($1 eq 'total');
        }
    }
    
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

print "</div>";

if (defined %states) {

    ## build the "image" div contents that represent the different states
    my $status_image = '';
    for my $status (sort keys %states) {
        ## don't include 'total' states
        next if ($status eq 'total');
    
        ## each status gives a percentage of the total command_count
        $status_image .= '<div class="status_bar_portion" style="width: ' . 
                           int( ($states{$status} / $command_count) * $progress_image_width) . 'px; background-color: ' . 
                           ($colors{$status} || 'rgb(0,0,0)') . ';">' . "</div>\n";
    }

    ## build the line that lists each status and its count
    my $status_list_line = '<li>states: ';

    for my $status (sort keys %states) {
        $status_list_line .= "$status (<span style='color:$colors{$status};'>$states{$status}</span>), ";
    }

    ## take off the trailing comma
    if ($status_list_line =~ /(.+)\,\s*$/) {
        $status_list_line = $1;
    }

    $status_list_line .= "</li>\n";

    print <<StateBox;
    <div id='summary'>
        <ul class='component'>
            <li><div class="component_progress_image">$status_image</div></li>
            <li>state: <span style='color: $colors{$overall_state}'>$overall_state</span></li>
            $status_list_line
       </ul>
    </div>
StateBox
}

print "</body></html>";

sub pageHeader {
    print <<heADerStuff;
<!DOCTYPE HTML PUBLIC "-//W3C//DTD XHTML 1.0 Strict//EN" --"http://www.w3.org/TR/xhtml1/DTD/xhtml1-strict.dtd"-->
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
        span.comment {
            background-color: rgb(200,200,200);
        }
        span.tag {
            color: rgb(100,100,100);
        }
        a {
            font-weight: bold;
        }
        li {
            list-style-type: none; 
            clear: left;
        }
        ul {
            width: 525px;
            padding-left: 5px;
            padding-bottom: 5px;
            margin-left: 20px;
        }
        ul.component {
            border: 1px solid grey;
            margin-bottom: 5px;
        }
        div.component_progress_image {
            border: 1px solid black;
            padding: 0px;
            margin: 0px;
            height: 10px;
            width: 500px;
            margin-top: 10px;
        }
        div.status_bar_portion {
            border: none;
            padding: 0px;
            margin: 0px;
            height: 10px;
            float: left;
        }
        #sourcecode {
            padding-left: 10px;
            margin-top: 10px;
            white-space: pre;
        }
    </style>
    <script type="text/javascript">
        function displaySummary() {
            // we want to get the summary and display it before the sourcecode div
            document.getElementById("sourcecode").parentNode.insertBefore(document.getElementById("summary"), document.getElementById("sourcecode"));
        }
    </script>
</head>
<body onLoad="javascript: displaySummary();">

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
