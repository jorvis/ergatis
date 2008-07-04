#!/usr/bin/perl -w

use strict;
use CGI;
use CGI::Carp qw(fatalsToBrowser);
use Ergatis::Common;
use Ergatis::SavedPipeline;
use HTML::Template;

my $q = new CGI;
print $q->header( -type => 'text/html' );

my $tmpl = HTML::Template->new( filename => 'templates/get_pipeline_templates.tmpl',
                                die_on_bad_params => 1,
                              );

my $project_templates = get_pipeline_templates( $q->param('path') );

$tmpl->param( PROJECT_TEMPLATES => $project_templates );

print $tmpl->output;

