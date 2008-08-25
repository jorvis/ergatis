#!/usr/bin/perl -w

use strict;
use CGI;
use CGI::Carp qw(fatalsToBrowser);
use Ergatis::Common;
use Ergatis::ConfigFile;
use Ergatis::SavedPipeline;
use HTML::Template;

my $q = new CGI;
print $q->header( -type => 'text/html' );

my $tmpl = HTML::Template->new( filename => 'templates/templates.tmpl',
                                die_on_bad_params => 1,
                              );

## read the ergatis config file
my $ergatis_cfg = new Ergatis::ConfigFile( -file => "ergatis.ini" );

my $global_templates = get_pipeline_templates( $ergatis_cfg->val('paths', 'global_saved_templates') );
my $recent_templates = get_pipeline_templates( $ergatis_cfg->val('paths', 'pipeline_build_area') );

my $projects = [];
## read each of the projects from the conf file
for my $project ( sort $ergatis_cfg->Parameters('projects') ) {
    my $project_data = {
        label               => $project,
        repository_root     => $ergatis_cfg->val('projects', $project ),
    };
    
    push @$projects, $project_data;
}

## 

$tmpl->param( GLOBAL_TEMPLATES    => $global_templates );
$tmpl->param( RECENT_TEMPLATES    => $recent_templates );
$tmpl->param( PROJECTS            => $projects );
$tmpl->param( QUICK_LINKS         => &get_quick_links($ergatis_cfg) );
$tmpl->param( SUBMENU_LINKS       => [
                                        #{ label => 'create project', is_last => 1, url => './create_project_form.cgi' },
                                     ] );

print $tmpl->output;

