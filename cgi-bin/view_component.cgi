#!/usr/local/bin/perl

use strict;

use XML::Twig;

use CGI qw(:standard);

use Tree::DAG_Node;
use Config::IniFiles;

my $conffile = param('conffile');

my $newcfg = new Config::IniFiles(-file => $conffile);
print header();	
print "<html><body><h3>$conffile</h3><br><pre>\n";
&print_config($newcfg);
print "</pre><br></body></html>";
exit;

sub print_config{
    my($cfg) = @_;
    my @sections = $cfg->Sections();
    foreach my $section (@sections){
	print "[$section]\n";
	my @parameters = $cfg->Parameters ($section);
	foreach my $param (@parameters){
	    my $value = $cfg->val($section,$param);
	    print "$param=$value\n";
	}
    }
}
