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
my $conf = param('conf');
my $location = param('location'); #first_child,before,after
my $name = param('name');


my $pipelineid = &get_pipeline_id($xmltemplate);
my $componentid = "$name"."_$$";

my $componentxml = XMLManip::get_component_xml($componentid,$pipelineid,$conf,$xmltemplate);
XMLManip::addxml($xmltemplate,$node,$componentxml,$location);

print redirect(-uri=>"./show_pipeline.cgi?xmltemplate=$xmltemplate&edit=1");

sub get_pipeline_id{
    my($file) = @_;
    my($id) = ($file =~ /pipeline(\d+)\./);
    if($id eq ""){
	$id = 0;
    }
    return $id;
}

