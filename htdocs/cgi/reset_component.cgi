#!/usr/bin/perl -w

use strict;
use CGI;
use File::Mirror;
use File::Basename;
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
        run_system_cmd("rm -rf $outputdir.old");
    }

    run_system_cmd("mv $outputdir $outputdir.old");
}

my $scratchdir =   $component_cfg->val('project', '$;TMP_DIR$;');
if(-d $scratchdir){
    if(-d "$scratchdir.old"){
        print STDERR "Removing $scratchdir.old\n";
        run_system_cmd("rm -rf $scratchdir.old");
    }

    run_system_cmd("mv $scratchdir $scratchdir.old");
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
                    # Set all serial/parallel nestings incomplete
                    # Valid to not have a state here so check if it exists before tyring to reset
                    if ( $elt->first_child('state') ) {
                        $elt->first_child('state')->set_text('incomplete');
                    }
                 }
             }
               },
             pretty_print => 'indented'
             );            
$twig->parse($pipeline_fh);
$twig->print_to_file("$pipeline.reset");
run_system_cmd("mv $pipeline.reset $pipeline");

## if we have a pipeline.xml.submitted file present we also want
## to get rid of it to
if (-e "$pipeline.submitted") {
    run_system_cmd("rm $pipeline.submitted");
}

if(-e "$component_xml"){
    run_system_cmd("mv $component_xml $component_xml.old");
}

## redirect to the template list page
print $q->redirect( -uri => url_dir_path($q) . "view_pipeline.cgi?instance=$pipeline" );
