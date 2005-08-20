#!/usr/local/bin/perl -w

use strict;
use CGI;
use CGI::Carp qw(fatalsToBrowser);
use Date::Manip;
use HTML::Template;
use XML::Twig;
use Data::Dumper;

my $q = new CGI;

print $q->header( -type => 'text/html' );

my $xml_input = $q->param("xml_input") || die "pass xml_input";

print_header();

if (-f $xml_input) {
    parseXMLFile($xml_input);
} else {
    die "don't know what to do with $xml_input";
}

print_footer();

exit(0);

sub parseXMLFile {
    my $file = shift;
    my $twig = new XML::Twig;

    $twig->parsefile($file);
   
    my $commandSetRoot = $twig->root;
    my $commandSet = $commandSetRoot->first_child('commandSet');
    parseCommandSet( $commandSet );
}


sub parseCommandSet {
    my ($commandSet) = @_;
   
    my $id = $commandSet->first_child('id')->text();
    
    print "<ul>\n";
    print "    <h1><span><b>commandSet</b> $id (" . $commandSet->first_child('configMapId')->text() . ")</span></h1>\n";
    
    my $cmds_found = 0;
    
    ## look at each child.
    for my $child ( $commandSet->children() ) {
        ## if it is a command, enter a li
        if ($child->gi eq 'command') {
            print "    <li>command " . $child->first_child('id')->text() . " (" . $child->first_child('configMapId')->text() . ")</li>\n";
            $cmds_found++;
            
        ## if it is a command set, parse through it.    
        } elsif ($child->gi eq 'commandSet') {
            parseCommandSet($child);
            $cmds_found++;
        }
    }
    
    ## if no commands or commandSets were found, check for a fileName reference,
    ##  which holds these in an external file.
    if (! $cmds_found) {
        my $fileName = $commandSet->first_child('fileName');
        
        if ($fileName) {
            if (-f $fileName->text) {
                parseXMLFile($fileName->text);
            } else {
#                print "<ul></ul>";
            }
        }
    }
    
    print "</ul>\n";
}


sub print_header {
    print <<HeAdER;

<html>

<head>
    <meta http-equiv="Content-Language" content="en-us">
    <meta http-equiv="Content-Type" content="text/html; charset=windows-1252">
    <title>pre-pre-alpha workflow interface prototype</title>
    <!-- <link rel="stylesheet" type="text/css" href="/concept/builder.css"> -->
    <style type="text/css">
        body {
            font-family: verdana, helvetica, arial, sans-serif;
            font-size: 10px;
            margin: 0;
            padding: 0;
        }
        
        ul h1 {
            font-weight: normal; 
            font-size: 100%; 
        }
        
        li {
            list-style-type: none; 
        }
        
    </style>
</head>

<body>

HeAdER
}

sub print_footer {
    print <<FooTER;
</body>

</html>
FooTER
}
