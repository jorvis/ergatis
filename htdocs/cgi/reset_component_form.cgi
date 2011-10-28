#!/usr/bin/perl -w

use strict;
use CGI;
use CGI::Carp qw(fatalsToBrowser);
use File::Basename;
use File::Spec;
use IO::File;
use Ergatis::Common;
use Ergatis::ConfigFile;

my $q = new CGI;

my $component = $q->param('component');
my $component_xml = $q->param('component_xml');
my $pipeline = $q->param('pipeline');
my $component_ini = $q->param('component_ini');

my $ergatis_cfg = new Ergatis::ConfigFile( -file => 'ergatis.ini' );
my $pipeline_url = $ENV{HTTP_REFERER};
my $cgi_script_dir = dirname( File::Spec->rel2abs( __FILE__ ) );
my $reset_script = "$cgi_script_dir/reset_component.cgi";

## Check whether or not we are logged-in and if so sudo to the current user 
## and execute our reset_component.cgi script as said user.
my $run_as = user_logged_in($ergatis_cfg);

## If we get an error we want to redirect to our built-in ergatis error page
my $error = { 'label' => 'last pipeline viewed', 'url' => $pipeline_url, 'cgi' => $q };

if ($run_as) {
    my $ergatis_cfg = new Ergatis::ConfigFile( -file => "ergatis.ini" );
    my $run_dir = $ergatis_cfg->val('paths', 'workflow_run_dir') || croak "workflow_run_dir not found in ergatis.ini";
    my $workflow_scripts_dir = "$run_dir/scripts";

    if (! -d $workflow_scripts_dir) {
        mkdir $workflow_scripts_dir || croak("Failed to create pipeline scripts directory $workflow_scripts_dir: $!");
    }

    chdir $run_dir || croak "Can't change to running directory $run_dir\n";

    my $reset_shell_script = "$workflow_scripts_dir/$run_as\_$component\_reset_" . time . ".sh";
    open (my $reset_fh, ">$reset_shell_script") || die "Cannot write reset component script $reset_shell_script: $!";

    print $reset_fh "#!/bin/bash", "\n\n";
    print $reset_fh "perl -I $cgi_script_dir/ " .
                    "$reset_script component=$component " .
                    "component_xml=$component_xml pipeline=$pipeline " .
                    "component_ini=$component_ini";
    close $reset_fh;                   
    chmod 0775, $reset_shell_script;

    run_system_cmd("sudo -u $run_as $reset_shell_script", $error);
} else{
    run_system_cmd("perl -I $cgi_script_dir/Ergatis $reset_script " .
                   "component=$component component_xml=$component_xml " .
                   "pipeline=$pipeline component_ini=$component_ini ",
  	  	   $error);
}

print $q->redirect( -uri => url_dir_path($q) . "view_pipeline.cgi?instance=$pipeline" );
