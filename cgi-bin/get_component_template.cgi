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
my $component_num;

if (defined $q->param('component_num')) {
    $component_num = $q->param('component_num');
} else {
    croak('need a component num');
}

my $shared_cfg = new Config::IniFiles( -file => "$repository_root/workflow_config_files/sharedconf.ini" );
my $workflowdocs_dir = $shared_cfg->val( 'init', '$;WORKFLOWDOCS_DIR$;' );

my $tmpl = HTML::Template->new( filename => 'templates/get_component_template.tmpl',
                                die_on_bad_params => 1,
                                global_vars => 1
                              );

my %selectable_labels = ( INPUT_FILE_LIST => 1, INPUT_FILE => 1, INPUT_DIRECTORY => 1 );

my $component_found = 0;
my $component_ini = "$workflowdocs_dir/${component_name}.config";
my $sections = [];

## make sure the configuration exists.
if ( -e $component_ini ) {
    
    ## read this config file
    my $component_cfg = new Config::IniFiles( -file => $component_ini );
    
    ## it's possible that later the first dd could be comments and the second the value
    for my $section ( $component_cfg->Sections() ) {
        ## config group is up to a whitepace
        $section =~ /^(.+?)\s+/;
        my $section_name = $1;

        ## any sections to skip?
        next if $section_name eq 'workflowdocs';
        next if $section_name eq 'include';
       
        push @$sections, { name => $section_name };
        
        
        for my $parameter ( $component_cfg->Parameters($section) ) {
            ## strip the parameter of the $; symbols
            my $label = $parameter;
            $label =~ s/\$\;//g;
            
            ## our variable naming 
            my $pretty_label = lc($label);
               $pretty_label =~ s/_/ /g;
            
            if ( exists $selectable_labels{$label} ) {
                push @{$$sections[-1]->{parameters}}, { label => $label, 
                                                        pretty_label => $pretty_label, 
                                                        value => $component_cfg->val($section, $parameter),
                                                        selectable => 1,
                                                        section => $section_name };            
            } else {
                push @{$$sections[-1]->{parameters}}, { label => $label, 
                                                        pretty_label => $pretty_label, 
                                                        value => $component_cfg->val($section, $parameter),
                                                        selectable => 0,
                                                        section => $section_name };
            }                
        }
    }
    
    $component_found++;
}

$tmpl->param( COMPONENT_FOUND => $component_found );
## $tmpl->param( COMPONENT_NAME    => $component_name );
$tmpl->param( SECTIONS => $sections );
$tmpl->param( COMPONENT_NUM => $component_num );

print $tmpl->output;

exit(0);
