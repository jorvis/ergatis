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
my $pipelinexml = param('pipelinexml');
my $location = param('location'); #first_child,before,after
my $id = param('id');


my $pipelineid = &get_pipeline_id($xmltemplate);

my $newpipelinexml = XMLManip::get_pipeline_xml($pipelineid,$xmltemplate,$pipelinexml);
XMLManip::addxml($xmltemplate,$node,$newpipelinexml,$location);

print redirect(-uri=>"./show_pipeline.cgi?xmltemplate=$xmltemplate&edit=1");

sub get_pipeline_id{
    my($file) = @_;
    my($id) = ($file =~ /pipeline(\d+)\./);
    if($id eq ""){
	$id = 0;
    }
    return $id;
}
