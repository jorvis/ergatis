#!/usr/local/bin/perl

use strict;

use XML::Twig;

use CGI qw(:standard);

use Tree::DAG_Node;
use Config::IniFiles;
use XMLManip;

my $root = param('root');

&create_directory($root);

my $xmltemplate = "$root/pipeline$$.xml";

my $newxml = XMLManip::get_root_commandset($xmltemplate);

my $t1 = new XML::Twig();
$t1->paste($newxml);

open FILE,"+>$xmltemplate" or die("Can't open file $xmltemplate");
$t1->print(\*FILE);
close FILE;

system("chmod 666 $xmltemplate");
system("chmod 666 $xmltemplate.ini");

print redirect(-uri=>"./show_pipeline.cgi?xmltemplate=$xmltemplate");

exit;

sub create_directory{
    my($dir) = @_;
    my $ret = system("mkdir -m 777 -p $dir");
    $ret >>= 8;
    return $ret;
}


