#!/usr/local/bin/perl -w

use strict;
use CGI;
use CGI::Carp qw(fatalsToBrowser);
use Config::IniFiles;
use HTML::Template;

my $q = new CGI;
print $q->header( -type => 'text/html' );

my $repository_root = $q->param('repository_root') || die "need a repository root";
my $component_name  = $q->param('component_name') || die "need a component name";
my $shared_cfg = new Config::IniFiles( -file => "$repository_root/workflow_config_files/sharedconf.ini" );
my $workflowdocs_dir = $shared_cfg->val( 'init', '$;WORKFLOWDOCS_DIR$;' );

my $tmpl = HTML::Template->new( filename => 'templates/get_component_template.tmpl',
                                die_on_bad_params => 1,
                              );

my $component_options = [];
my $component_found = 0;

my $component_ini = "$workflowdocs_dir/${component_name}conf.ini";

## make sure the configuration exists.
if ( -e $component_ini ) {
    
    ## read this config file
    my $component_cfg = new Config::IniFiles( -file => $component_ini );

    ## it's possible that later the first dd could be comments and the second the value
    for my $section ( $component_cfg->Sections() ) {
        for my $parameter ( $component_cfg->Parameters($section) ) {
            ## strip the parameter of the $; symbols
            my $label = $parameter;
            $label =~ s/\$\;//g;

            push @{$component_options}, { label => $label, value => ($component_cfg->val($section, $parameter) || '') };
        }
    }
    
    $component_found++;
}

$tmpl->param( COMPONENT_FOUND   => $component_found );
$tmpl->param( COMPONENT_NAME    => $component_name );
$tmpl->param( COMPONENT_OPTIONS => $component_options );

print $tmpl->output;

exit(0);
