#!/usr/bin/perl -w

use strict;
use CGI;
use CGI::Carp qw(fatalsToBrowser);
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
       	#print STDERR `rm -rf $outputdir.old`;
    }

    run_system_cmd("mv $outputdir $outputdir.old");
    #print STDERR `mv $outputdir $outputdir.old`;
}

my $scratchdir =   $component_cfg->val('project', '$;TMP_DIR$;');
if(-d $scratchdir){
    if(-d "$scratchdir.old"){
        print STDERR "Removing $scratchdir.old\n";
        run_system_cmd("rm -rf $scratchdir.old");
        #print STDERR `rm -rf $scratchdir.old`;
    }

    run_system_cmd("mv $scratchdir $scratchdir.old");
    #print STDERR `mv $scratchdir $scratchdir.old`;
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
run_system_cmd("mv $pipeline.reset $pipeline");
#print STDERR `mv $pipeline.reset $pipeline`;

if(-e "$component_xml"){
    run_system_cmd("mv $component_xml $component_xml.old");
    #print STDERR `mv $component_xml $component_xml.old`;
}

## If any iterator directories exist we want to rename them to <iterator directory>.old 
## to make sure that the iterator XML's from the previous run don't show up in the re-run.
my @iterator_dirs = get_iterator_directories( dirname($component_xml) );
foreach my $iter_dir (@iterator_dirs) {
    run_system_cmd("mv $iter_dir $iter_dir.old");
    #print STDERR `mv $iter_dir $iter_dir.old`;
}

## redirect to the template list page
print $q->redirect( -uri => url_dir_path($q) . "view_pipeline.cgi?instance=$pipeline" );

#############################################################
#                     SUBROUTINES                           #
#############################################################

###
# Runs a system command wrapped in an eval to catch any issues that 
# could occur during the process. If an error occurs this script will die
# spitting out an error.
###
sub run_system_cmd {
    my $cmd = shift;
    my @cmd_output;

    eval {
        @cmd_output = qx{$cmd 2>&1};
        if ( ($? << 8) != 0 ) {
            die "@cmd_output";
        }
    }; 
    if ($@) {
        die "Error executing command $cmd: $@";
    }
}

###
# Gets all iterator directories in the component output directory
###
sub get_iterator_directories {
    my $component_dir = shift;

    opendir (DIR, $component_dir) || die "Cannot scan directory $component_dir: $!";
    my @iter_dirs = grep { /i\d+$/  && -d "$component_dir/$_" } readdir(DIR);
    @iter_dirs = map { $component_dir . "/" . $_ } @iter_dirs;
    closedir DIR;
            
    return @iter_dirs;     
}
