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

my $ergatis_cfg = new Ergatis::ConfigFile( -file => "ergatis.ini" );
my $auth_method = $ergatis_cfg->val('authentication', 'authentication_method');
my $username = user_logged_in($ergatis_cfg);

my $instance = $q->param('instance');
my $repository_root = $q->param('repository_root');
my $template_name = "temp$$";
my ($pipeline_dir) = ($instance =~ /(.*\/pipeline\/[A-Z0-9]+\/)pipeline.xml/);
my ($pipeline_id) = ($instance =~ /.*\/pipeline\/([A-Z0-9]+)\//);
my $pipeline_layout = $pipeline_dir."pipeline.layout";
my @pipeline_configs = glob($pipeline_dir."*.config");

unless ($auth_method eq 'open' || defined($username)) {
    print $q->header( -type => 'text/html' );
    print_error_page( ergatis_cfg => $ergatis_cfg,
                      message => "You must be logged in to clone pipelines",
                      links => [ 
                                 { label => "pipeline list", is_last => 1, url => "./pipeline_list.cgi?repository_root=$repository_root" },
                               ],
                    );
    exit(0);
}

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
