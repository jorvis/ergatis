#!/usr/local/bin/perl -w

use strict;
use CGI;
use CGI::Carp qw(fatalsToBrowser);
use Config::IniFiles;
use HTML::Template;

my $q = new CGI;
print $q->header( -type => 'text/html' );

my $tmpl = HTML::Template->new( filename => 'templates/index.tmpl',
                                die_on_bad_params => 1,
                              );

## read the ergatis config file
my $ergatis_cfg = new Config::IniFiles( -file => "ergatis.ini" );


## build the project list
my $registered_projects = [];
for my $label ( sort $ergatis_cfg->Parameters('projects') ) {
    push @$registered_projects, { 
                                    label => $label,
                                    repository_root => $ergatis_cfg->val('projects', $label),
                                };
}

$tmpl->param( REGISTERED_PROJECTS => $registered_projects );

print $tmpl->output;
