#!/usr/local/bin/perl

use strict;

use XML::Twig;

use CGI qw(:standard);

use Tree::DAG_Node;

#edit fields: name,configMapId,state,type

my $xmltemplate = param('xmltemplate');
my $node = param('node');

my $dowrite = param('dowrite');
my $editname = param('editname');
my $editconfigmapid = param('editconfigmapid');
my $editstate = param('editstate');
my $edittype = param('edittype');

if($dowrite == 1){
    my $t1 = new XML::Twig( TwigHandlers => { 'commandSet' =>
					      sub {
						  my ($t, $elt) = @_;
						  my $configmapid = $elt->first_child("configMapId");
						  my $id = $elt->att("id");

						  if(((defined $configmapid) && ($configmapid->text() eq $node)) || ($id eq $node)){
						      &set_node($elt,"name",$editname);
						      &set_node($elt,"configMapId",$editconfigmapid);
						      &set_node($elt,"state",$editstate);
						      &set_node($elt,"type",$edittype);
						  }
					      }
				      },
			pretty_print => 'indented'
			);

    $t1->parsefile($xmltemplate);			
    open FILE, "+>$xmltemplate" or die("Can't open output xml file $xmltemplate.new : $!");
    $t1->print(\*FILE);
    close FILE;
    `chmod 666 $xmltemplate`;

    print redirect(-uri=>"./show_pipeline.cgi?xmltemplate=$xmltemplate&node=$node");
}
else{
    print header();
    print "<html><body><a href='show_pipeline.cgi?xmltemplate=$xmltemplate'>[view workflow]</a>&nbsp;<a href='remove_commandset.cgi?xmltemplate=$xmltemplate&node=$node'>[delete]</a><br>\n";
    print "<form><input type='hidden' name='xmltemplate' value='$xmltemplate'><input type='hidden' name='node' value='$node'><input type='hidden' name='dowrite' value='1'>\n";
    my $t1 = new XML::Twig( TwigHandlers => { 'commandSet' =>
						  sub {
						      my ($t, $elt) = @_;
						      my $configmapid = $elt->first_child("configMapId");
						      my $id = $elt->att("type");
						      
						      my $newxml;
						      if(((defined $configmapid) && ($configmapid->text() eq $node)) || ($id eq $node)){
							  print "Name: <input type='text' name='editname' value='",$elt->first_child_text("name"),"'><br>\n";
							  print "configMapId: <input type='text' name='editconfigmapid' value='",$elt->first_child_text("configMapId"),"'><br>\n";
							  print "state: <input type='text' name='editstate' value='",$elt->first_child_text("state"),"'><br>\n";
							  print "type: <input type='text' name='edittype' value='",$elt->first_child_text("type"),"'><br>\n";
						      }
						  }
					      });
    $t1->parsefile($xmltemplate);
    print "<input type='submit' value='submit'></body></html>\n";
}



sub set_node{
    my($elt,$eltname,$editvalue) = @_;
    if($editvalue ne ""){
	if($elt->first_child($eltname)){
	    $elt->first_child($eltname)->set_text($editvalue);
	}
	else{
	    my $newxml = parse XML::Twig::Elt("<$eltname>$editvalue</$eltname>");
	    $newxml->paste('last_child',$elt);
	}
    }
}
    
