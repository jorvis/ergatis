#!/usr/bin/perl -w

use strict;
use CGI;
use CGI::Carp qw(fatalsToBrowser);
use Ergatis::Common;
use Ergatis::ConfigFile;
use File::Mirror;
use HTML::Template;
use File::Basename;

my $q = new CGI;

my $instance = $q->param('instance');
my $repository_root = $q->param('repository_root');
my $template_name = "temp$$";
my ($pipeline_dir) = ($instance =~ /(.*\/pipeline\/\d+\/)pipeline.xml/);
my ($pipeline_id) = ($instance =~ /.*\/pipeline\/(\d+)\//);
my $pipeline_layout = $pipeline_dir."pipeline.layout";
my @pipeline_configs = glob($pipeline_dir."*.config");


## if we get to here all needed options were passed
mkdir('/tmp/pipelines_building');
my $output_template = "/tmp/pipelines_building/cloneof_$pipeline_id";
 
$output_template =~ s/\s/_/g;
 
## copy the template
mkdir($output_template);
mirror($pipeline_layout, "$output_template/pipeline.layout" );
foreach my $config (@pipeline_configs){
    my $configfile = basename($config,"*.config");
    print STDERR "CONFIG $config $output_template/$configfile\n";

    mirror($config,"$output_template/$configfile");
}

## redirect to the template list page
print $q->redirect( -uri => url_dir_path($q) . "build_pipeline.cgi?autoload_template=$output_template&repository_root=".$q->param('repository_root') );
