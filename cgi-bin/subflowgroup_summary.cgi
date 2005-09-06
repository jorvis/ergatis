#!/usr/local/bin/perl -w

use strict;
use CGI;
use CGI::Carp qw(fatalsToBrowser);
use File::Basename;
use XML::Twig;

my $q = new CGI;

print $q->header( -type => 'text/html' );

my $xml_input = $q->param("xml_input") || die "pass xml_input";

## the subflowNgroupsN.xml file is always small, so we can DOM parse
my $sg_twig = XML::Twig->new();
   $sg_twig->parsefile($xml_input);

my $parent_commandset = $sg_twig->root->first_child('commandSet');

## the subflowNgroupsN.xml file will only have one command
my $command = $parent_commandset->first_child('command');

print "    <h1>analysis group elements</h1>\n";

## we need to parse through the params to get the output xml
##  perhaps later we can build the command as well
## the subflow_subflowNgroupsN.xml file
my $commandset = $parent_commandset->first_child('commandSet');
my $ssg_file = $commandset->first_child('fileName')->text();

my $ssg_twig = XML::Twig->new( twig_roots => {
                                    'commandSet/commandSet' => \&process_iterated_commandset,
                               },
               );
$ssg_twig->parsefile($ssg_file);



sub process_iterated_commandset {
    my ($twig, $commandset) = @_;
    
    my $file = $commandset->first_child('fileName')->text();
    my $name = basename($file, '.xml');
    my $state = $commandset->first_child('state')->text();
    
    print <<SubFlowFile;
    <div id='${name}_bar' class='subflowbar'>
        <div class='leftside'>
            <img id='${name}_arrow' class='expander' src='/cram/arrow_right.gif' onclick='toggle_subflow_display("$name", "$file");' alt='expand' title='expand'>
            <img id='${name}_img' class='status' src='/cram/status_$state.png' title='$state' alt='$state'>
            $name
        </div>
        <div class='rightside'>
            <img class='reloader' src='/cram/reload_blue.png' onclick='reload_subflow("$name", "$file")' alt='reload' title='reload'>
        </div>
    </div>
    <div id='${name}_data' class='subflowdata' style='display: none;'></div>
SubFlowFile

    $twig->purge;
}
