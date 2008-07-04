#!/usr/bin/perl -w

use strict;
use CGI;
use CGI::Carp qw(fatalsToBrowser);
use Ergatis::Common;
use Ergatis::ConfigFile;
use File::Mirror;
use HTML::Template;

my $q = new CGI;

my $input_template = $q->param('input_template');
my $repository_root = $q->param('repository_root');
my $template_name = $q->param('template_name');

my $ergatis_cfg = new Ergatis::ConfigFile( -file => "ergatis.ini" );

unless ( $input_template && $repository_root && $template_name ) {
    print $q->header( -type => 'text/html' );
    
    print_error_page( ergatis_cfg => $ergatis_cfg,
          message => "input_template, repository_root and template_name are required options to this script. " .
                     "This should not happen, please file a bug report at ergatis.sf.net.",
          links => [ 
                        { label => "bug report page", 
                          is_last => 1, 
                          url => "http://sf.net/projects/ergatis" },
                   ],
    );
    
    exit();
}

## if we get to here all needed options were passed
my $output_template = "$repository_root/workflow/project_saved_templates/$template_name";
 
$output_template =~ s/\s/_/g;
 
## copy the template
mirror( $input_template, $output_template );

## redirect to the template list page
print $q->redirect( -uri => url_dir_path($q) . "templates.cgi" );
