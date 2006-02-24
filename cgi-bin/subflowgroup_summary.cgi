#!/usr/local/bin/perl -w

use strict;
use CGI;
use CGI::Carp qw(fatalsToBrowser);
use File::Basename;
use XML::Twig;

my $q = new CGI;

print $q->header( -type => 'text/html' );

my $xml_input = $q->param("xml_input") || die "pass xml_input";

my $xml_input_fh;
if ($xml_input =~ /\.gz/) {
    open($xml_input_fh, "<:gzip", "$xml_input") || die "can't read $xml_input: $!"; 
} elsif ( ! -e $xml_input && -e "$xml_input.gz" ) {
    open($xml_input_fh, "<:gzip", "$xml_input.gz") || die "can't read $xml_input: $!"; 
} else {
    open($xml_input_fh, "<$xml_input") || die "can't read $xml_input: $!";       
}

## the subflowNgroupsN.xml file is always small, so we can DOM parse
my $sg_twig = XML::Twig->new();
   $sg_twig->parse($xml_input_fh);

my $parent_commandset = $sg_twig->root->first_child('commandSet');

## the subflowNgroupsN.xml file will only have one command
my $command = $parent_commandset->first_child('command');

print "    <h1>analysis group elements</h1>\n";

## we need to parse through the params to get the output xml
##  perhaps later we can build the command too
## the subflow_subflowNgroupsN.xml file
my $commandset = $parent_commandset->first_child('commandSet');
my $ssg_file = $commandset->first_child('fileName')->text();

my $ssg_file_fh;
if ($ssg_file =~ /\.gz/) {
    open($ssg_file_fh, "<:gzip", "$ssg_file") || die "can't read $ssg_file: $!"; 
} elsif ( ! -e $ssg_file && -e "$ssg_file.gz" ) {
    open($ssg_file_fh, "<:gzip", "$ssg_file.gz") || die "can't read $ssg_file: $!"; 
} else {
    open($ssg_file_fh, "<$ssg_file") || die "can't read $ssg_file: $!";       
}

my $ssg_twig = XML::Twig->new( twig_roots => {
                                    'commandSet/commandSet' => \&process_iterated_commandset,
                               },
               );
$ssg_twig->parse($ssg_file_fh);



sub process_iterated_commandset {
    my ($twig, $commandset) = @_;
    
    my $file = $commandset->first_child('fileName')->text();
    my $name = basename($file, ('.xml', '.gz'));
    my $state = $commandset->first_child('state')->text();
    
    ## check for any messages
    my $msg_html = '';
    if ( $commandset->has_child('status') ) {
        if ( $commandset->first_child('status')->has_child('message') ) {
            my $msg = $commandset->first_child('status')->has_child('message')->text();
            
            ## don't display it if it is just a 'finished' message
            unless ( $msg =~ /Command set with name.+ finished/i ) {
                $msg_html = "<div class='messageblock'>$msg</div>\n";
            }
        }
    }
    
    print <<SubFlowFile;
    <div id='${name}_bar' class='subflowbar'>
        <div class='leftside'>
            <img id='${name}_arrow' class='expander' src='/ergatis/arrow_right.gif' onclick='toggle_subflow_display("$name", "$file");' alt='expand' title='expand'>
            <img id='${name}_img' class='status' src='/ergatis/status_$state.png' title='$state' alt='$state'>
            $name
        </div>
        <div class='rightside'>
            <img class='reloader' src='/ergatis/reload_blue.png' onclick='reload_subflow("$name", "$file")' alt='reload' title='reload'>
        </div>
    </div>
    $msg_html
    <div id='${name}_data' class='subflowdata' style='display: none;'></div>
SubFlowFile

    $twig->purge;
}
