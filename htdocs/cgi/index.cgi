#!/usr/bin/perl -w

use strict;
use CGI;
use CGI::Carp qw(fatalsToBrowser);
use Data::Dumper;
use Ergatis::Common;
use Ergatis::ConfigFile;
use Ergatis::Monitor;
use HTML::Template;
use Storable;
use XML::Twig;

my $q = new CGI;
print $q->header( -type => 'text/html' );

umask(0000);

## this toggle will force a rescan of the pipelines rather than pulling from 
##  the storable object.
my $update_cache = $q->param('update_cache') || 0;

my $tmpl = HTML::Template->new( filename => 'templates/index.tmpl',
                                die_on_bad_params => 0,
                              );

## read the ergatis config file
my $ergatis_cfg = new Ergatis::ConfigFile( -file => "ergatis.ini" );

## build the project list
my $registered_projects = [];
for my $label ( sort $ergatis_cfg->Parameters('projects') ) {
    push @$registered_projects, { 
                                    label => $label,
                                    repository_root => $ergatis_cfg->val('projects', $label),
                                };
}

$tmpl->param( REGISTERED_PROJECTS => $registered_projects );
$tmpl->param( UPDATE_CACHE        => $update_cache );
$tmpl->param( DEFAULT_PROJECT_ROOT => $ergatis_cfg->val( 'paths', 'default_project_root') || '' );
$tmpl->param( QUICK_LINKS         => &get_quick_links($ergatis_cfg) );
$tmpl->param( SUBMENU_LINKS       => [
                                        { label => 'update cache', is_last => 1, url => './index.cgi?update_cache=1' },
                                     ] );

print $tmpl->output;

