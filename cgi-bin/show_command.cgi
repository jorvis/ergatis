#!/usr/local/bin/perl

use strict;

use XML::Twig;

use CGI qw(:standard);

use Tree::DAG_Node;

#edit fields: name,configMapId,state,type

my $xmltemplate = param('xmltemplate');
my $node = param('node');

print header();

print "<html><body><a href='show_pipeline.cgi?xmltemplate=$xmltemplate'>[view sub-workflow]</a>&nbsp;<a href='javascript:back();'>[back]</a>";

my $t1 = new XML::Twig( TwigHandlers => { 'command' =>
					      sub {
						  my ($t, $elt) = @_;
						  my $configmapid = $elt->first_child("configMapId");
						  my $id = $elt->first_child_text("id");

						  my $runstr = &runstring($elt);

						 

						  if(((defined $configmapid) && ($configmapid->text() eq $node)) || ($id eq $node)){
						      print "<br><pre>$runstr</pre><br><hr>";
						      print "<pre>";
						      my $text = $elt->sprint();
						      my @lines  = split(/\n/,$text);
						      foreach my $line (@lines){
							  if($line =~ /\/.*\.xml/){
							      $line =~ s/(\/.*\.xml)/<a href='show_pipeline.cgi?&xmltemplate=$1'>$1<\/a>/mg;
							  }
							  else{
							      $line =~ s/(\/.*\.\w+)/<a href='show_file.cgi?&file=$1'>$1<\/a>/mg;
							  }
							  print $line,"\n";
						      }
						      print "</pre>";
						  }
					      }
				      },
			output_filter => 'html',
			pretty_print => 'indented');

$t1->parsefile($xmltemplate);

print "</body></html>";
	
sub runstring{
    my($elt) = @_;
    my @params = $elt->children('param');
    my @args = $elt->children('arg');
   
    my @allparams;
    my @allargs;
    
    my $command;
    my $stderr;
    my $stdout;

    foreach my $param (@params){
	my $keystr = $param->first_child('key')->text();
	my $valuestr = $param->first_child('value')->text();
	if($keystr eq 'command'){
	    $command = $valuestr;
	}
	elsif($keystr eq 'stderr'){
	    $stderr = $valuestr;
	}
	elsif($keystr eq 'stdout'){
	    $stdout = $valuestr;
	}
	elsif($keystr =~ /^-/){
	    push @allparams, "$keystr=$valuestr";
	}
	else{
	    push @allparams, "--$keystr=$valuestr";
	}
    }
    foreach my $arg (@args){
	my $argstr = $arg->text();
	push @allargs,$argstr;
    }

    my $execstr = "$command ".join(" ",@allparams)." ".join(" ",@allargs);
    if($stdout ne ""){
	$execstr = "($execstr";
	$execstr .= " > $stdout)";
    }
    if($stderr ne ""){
	$execstr .= " >& $stderr";
    }
    return $execstr;
}
    
    
