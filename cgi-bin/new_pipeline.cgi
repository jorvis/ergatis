#!/usr/local/bin/perl

use strict;

use XML::Twig;

use CGI qw(:standard);

use Tree::DAG_Node;
use Config::IniFiles;

my $root = param('root');

umask 000;

## create the directory for the pipeline
if (create_directory("$root/$$")) {
    die "couldn't create directory $root/$$ for pipeline\n";
}

my $xmltemplate = "$root/$$/pipeline.xml";

my $newxml = &get_skeleton_commandset();

open FILE,"+>$xmltemplate" or die("Can't open file $xmltemplate");
print FILE $newxml;
close FILE;

open FILE,"+>$xmltemplate.ini" or die("Can't open file $xmltemplate");
print FILE "[start]\n\n";
close FILE;

system("chmod 666 $xmltemplate");

print redirect(-uri=>"./show_pipeline.cgi?xmltemplate=$xmltemplate");

exit;

sub get_skeleton_commandset{
    return '<commandSetRoot xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:schemaLocation="commandSet.xsd">
            <commandSet type="serial">
                <configMapId>start</configMapId>
            </commandSet>
            </commandSetRoot>';
}


sub create_directory{
    my($dir) = @_;
    my $ret = system("mkdir -m 777 -p $dir");
    $ret >>= 8;
    return $ret;
}


