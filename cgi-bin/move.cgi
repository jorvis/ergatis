#!/usr/local/bin/perl

use strict;

use XML::Twig;

use CGI qw(:standard);
#use CGI::Carp qw(fatalsToBrowser set_message);

use Tree::DAG_Node;
use XMLManip;

my $xmltemplate = param('xmltemplate');
my $node = param('node');
my $location = param('location');
my $target = param('target');

my $xml = XMLManip::movexml($xmltemplate,$node,$target,$location);


print redirect(-uri=>"./show_pipeline.cgi?xmltemplate=$xmltemplate&edit=1");
