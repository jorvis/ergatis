#!/usr/local/bin/perl

use strict;

use XML::Twig;

use CGI qw(:standard);

use Tree::DAG_Node;

#edit fields: name,configMapId,state,type

my $xmltemplate = param('xmltemplate');
my $node = param('node');

print header();

print "<html><body><a href='show_pipeline.cgi?xmltemplate=$xmltemplate'>[view workflow]</a>";

my $t1 = new XML::Twig( TwigHandlers => { 'commandSet' =>
					      sub {
						  my ($t, $elt) = @_;
						  my $configmapid = $elt->first_child("configMapId");
						  my $id = $elt->first_child_text("id");
						  my $parentxml = $elt->first_child_text("parentFileName");

						  if($parentxml){
						      print "<a href='show_pipeline.cgi?xmltemplate=$parentxml'>[view parent workflow]</a><br>";
						  }
						  if(((defined $configmapid) && ($configmapid->text() eq $node)) || ($id eq $node)){
						      print "<pre>";
						      $elt->print();
						      print "</pre>";
						  }
					      }
				      },
			output_filter => 'html',
			pretty_print => 'indented');

$t1->parsefile($xmltemplate);

print "</body></html>";
	
