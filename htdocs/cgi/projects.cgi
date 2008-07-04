#!/usr/bin/perl -w

use strict;
use CGI;
use CGI::Carp qw(fatalsToBrowser);
use Ergatis::Common;
use Ergatis::ConfigFile;
use HTML::Template;

my $q = new CGI;
print $q->header( -type => 'text/html' );

my $tmpl = HTML::Template->new( filename => 'templates/projects.tmpl',
                                die_on_bad_params => 0,
                              );

## read the ergatis config file
my $ergatis_cfg = new Ergatis::ConfigFile( -file => "ergatis.ini" );

my $projects = [];
## read each of the projects from the conf file
for my $project ( $ergatis_cfg->Parameters('projects') ) {
    my $project_data = {
        label               => $project,
        project_url         => "./pipeline_list.cgi?repository_root=" . $ergatis_cfg->val('projects', $project ),
        repository_root     => $ergatis_cfg->val('projects', $project ) || '?',
        ergatis_dir         => '?',
    };
    
    my $shared_cfg_path = "$$project_data{repository_root}/workflow/project.config";
    
    ## make sure it exists
    if ( -e $shared_cfg_path ) {
        my $shared_cfg = new Ergatis::ConfigFile( -file => $shared_cfg_path );

        if ( defined $shared_cfg ) {
            $$project_data{ergatis_dir} = $shared_cfg->val('project', '$;ERGATIS_DIR$;');
        }
    }

    push @$projects, $project_data;
}

## 

$tmpl->param( PROJECTS            => $projects );
$tmpl->param( QUICK_LINKS         => &get_quick_links($ergatis_cfg) );
$tmpl->param( SUBMENU_LINKS       => [
                                        #{ label => 'create project', is_last => 1, url => './create_project_form.cgi' },
                                     ] );

print $tmpl->output;

