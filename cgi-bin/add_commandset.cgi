#!/usr/local/bin/perl

use strict;

use XML::Twig;

use CGI qw(:standard);
#use CGI::Carp qw(fatalsToBrowser set_message);

use Tree::DAG_Node;
use Config::IniFiles;
use XMLManip;

my $xmltemplate = param('xmltemplate');
my $node = param('node');
my $type = param('type');
my $location = param('location'); #first_child,before,after
my $name = param('name');

my $csid = "$$";

my $csxml = XMLManip::get_commandset_xml($type,$csid,$xmltemplate);
XMLManip::addxml($xmltemplate,$node,$csxml,$location);

print redirect(-uri=>"./show_pipeline.cgi?xmltemplate=$xmltemplate&edit=1");

