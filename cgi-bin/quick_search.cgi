#!/usr/local/bin/perl -w

use strict;
use CGI;
use CGI::Carp qw(fatalsToBrowser);
use Ergatis::Common;
use Ergatis::ConfigFile;
use HTML::Template;

my $q = new CGI;

## this template will only be used if we can't find anything.
my $tmpl = HTML::Template->new( filename => 'templates/quick_search.tmpl',
                                die_on_bad_params => 1,
                              );
my $crit = $q->param('quick_search_crit') || '';

## read the ergatis config file
my $ergatis_cfg = new Ergatis::ConfigFile( -file => "ergatis.ini" );

## search projects:*, paths:default_project_root, and pipeline IDs for each if it's numeric
for my $project ( $ergatis_cfg->Parameters('projects') ) {
    if ( lc($project) eq lc($crit) ) {
        ## redirect to project pipeline list and stop
        print $q->redirect( -uri => url_dir_path($q) . "pipeline_list.cgi?repository_root=" . $ergatis_cfg->val('projects', $project) );
        exit(0);
    }
    
    ## if it's numeric search the pipelines under this project
    if ( $crit =~ /^\d+$/ ) {
        ## search pipelines
        my $pipeline_dir = $ergatis_cfg->val('projects', $project) . '/workflow/runtime/pipeline';
        if (-d $pipeline_dir) {
            opendir(my $pdh, $pipeline_dir) || die "failed to open pipeline directory $pipeline_dir";
            while (my $dir = readdir $pdh) {
                if ($dir eq $crit && -e "$pipeline_dir/$dir/pipeline.xml") {
                    print $q->redirect( -uri => url_dir_path($q) . "view_pipeline.cgi?instance=$pipeline_dir/$dir/pipeline.xml" );
                    exit(0);
                }
            }
        }
    }
}

print $q->header( -type => 'text/html' );
$tmpl->param( QUICK_LINKS         => &get_quick_links($ergatis_cfg) );
$tmpl->param( SUBMENU_LINKS       => [
                                        #{ label => 'create project', is_last => 1, url => './create_project_form.cgi' },
                                     ] );

print $tmpl->output;

