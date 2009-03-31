#!/usr/bin/perl -w

use strict;
use CGI;
use CGI::Carp qw(fatalsToBrowser);
use File::Mirror;
use XML::Twig;
use Ergatis::Common;
use Ergatis::ConfigFile;
use Ergatis::Pipeline;


my $q = new CGI;

my $component = $q->param('component');
my $component_xml = $q->param('component_xml');
my $pipeline = $q->param('pipeline');
my $component_ini = $q->param('component_ini');

my $pipeline_fh;
if ($pipeline =~ /\.gz/) {
    open($pipeline_fh, "<:gzip", "$pipeline") || die "can't read $pipeline: $!"; 
} elsif ( -e "$pipeline.gz" ) {
    open($pipeline_fh, "<:gzip", "$pipeline.gz") || die "can't read $pipeline: $!"; 
} else {
    open($pipeline_fh, "<$pipeline") || die "can't read $pipeline: $!";       
}


my $component_cfg = new Ergatis::ConfigFile( -file => $component_ini );

#clear output and scratch dirs
my $outputdir =   $component_cfg->val('output', '$;OUTPUT_DIRECTORY$;');
if(-d $outputdir){
    if(-d "$outputdir.old"){
	print STDERR "Removing $outputdir.old\n";
       	print STDERR `rm -rf $outputdir.old`;
    }
    print STDERR `mv $outputdir $outputdir.old`;
}

my $scratchdir =   $component_cfg->val('project', '$;TMP_DIR$;');
if(-d $scratchdir){
    if(-d "$scratchdir.old"){
	print STDERR "Removing $scratchdir.old\n";
	print STDERR `rm -rf $scratchdir.old`;
    }
    print STDERR `mv $scratchdir $scratchdir.old`;
}

my $twig = new XML::Twig(
			 twig_handlers =>         
			 { commandSet => sub {
			     my( $twig, $elt)= @_;
			     my $name  = $elt->first_child( 'name')->text;
			     if($name eq $component){
				 #Set init commands for the component incomplete
				 my @commands = $elt->children('command');
				 foreach my $command (@commands){
				     if($command->first_child('state')){
					 $command->first_child('state')->set_text('incomplete');
				     }
				 }
				 #Set the rest of the component incomplete
				 $elt->first_child('state')->set_text('incomplete');
				 $elt->first_child('commandSet')->first_child('state')->set_text('incomplete');
			     }
			     elsif($name eq 'serial' || $name eq 'parallel' || $name eq "start pipeline:"){
				 #Set all serial/parallel nestings incomplete
				 $elt->first_child('state')->set_text('incomplete');
			     }
			 }
		       },
			 pretty_print => 'indented'
			 );            
$twig->parse($pipeline_fh);
$twig->print_to_file("$pipeline.reset");
print STDERR `mv $pipeline.reset $pipeline`;
if(-e "$component_xml"){
    print STDERR `mv $component_xml $component_xml.old`;
}

## redirect to the template list page
print $q->redirect( -uri => url_dir_path($q) . "view_pipeline.cgi?instance=$pipeline" );
