#!/usr/bin/perl -w

=head1 DESCRIPTION

There are three use cases for component config retrieval.

Arguments needed depend on the mode used.

1. completely new component instance
    - mode: new
    - default configuration
    - not automatically saved

2. clone of a component while building a pipeline
    - mode: clone
    - exact copy of component's values except OUTPUT_TOKEN
    - not automatically saved

3. copy of component while pipeline is copied or imported
    - mode: exact
    - exact copy of component, including OUTPUT_TOKEN
    - automatically saved

=cut

use strict;
use CGI;
use CGI::Carp qw(fatalsToBrowser);
use Ergatis::ConfigFile;
use HTML::Template;
use File::Basename;

my $q = new CGI;
print $q->header( -type => 'text/html' );

my $mode = $q->param('mode') || croak("need a mode");
my $component_ini  = $q->param('path') || croak("need a path");
my $component_id    = $q->param('component_id') || croak("need a component id");

my $tmpl = HTML::Template->new( filename => 'templates/get_component_template.tmpl',
                                die_on_bad_params => 1,
                                global_vars => 1
                              );

my %selectable_labels = ( INPUT_FILE_LIST => 1, INPUT_FILE => 1, INPUT_DIRECTORY => 1,
                          QUERY_BSML_FILE_LIST => 1, QUERY_BSML_FILE => 1, QUERY_BSML_DIRECTORY => 1,
                          MATCH_BSML_FILE_LIST => 1, MATCH_BSML_FILE => 1, MATCH_BSML_DIRECTORY => 1,
                          JACCARD_OUTPUT_LIST => 1, BLAST_STORED_FILE => 1, DB_FILTER_LIST => 1,
                          PANGENOME_INPUT_TABLE => 1,);

my $component_found = 0;
my %sections = ( basic => [], advanced => [] );

## make sure the configuration exists.
if ( -e $component_ini ) {
    
    ## read this config file
    my $component_cfg = new Ergatis::ConfigFile( -file => $component_ini );

    ## A configuration file that contains a set of common parameters that should
    ## be in every component
    my $common_ini = dirname($component_ini) . "/shared_parameters.config";
    if (-e $common_ini) {
        my $common_cfg = new Ergatis::ConfigFile( -file => $common_ini );
        $component_cfg->merge_configs($common_cfg);
    }

    ## it's possible that later the first dd could be comments and the second the value
    for my $section ( $component_cfg->Sections() ) {
        my $section_type = 'basic';
        if ($section =~ /workflowdocs|include|component|interface|dce/) {
            $section_type = 'advanced';
        }
       
        push @{$sections{$section_type}}, { name => $section };
        
        for my $parameter ( $component_cfg->Parameters($section) ) {
            ## strip the parameter of the $; symbols
            my $label = $parameter;
            $label =~ s/\$\;//g;
            
            ## our variable naming 
            my $pretty_label = lc($label);
               $pretty_label =~ s/_/ /g;
            
            my $value;
            
            ## when cloning, OUTPUT_TOKEN isn't copied
            if ( $mode eq 'clone' && $parameter eq '$;OUTPUT_TOKEN$;' ) {
                $value = 'default';
            } else {
                $value = $component_cfg->val($section, $parameter);
            }
            my @enumeration;
	    if($value =~ /^\s*\[(.*)\]/){
		my @vals = split(/\|/,$1);
		foreach my $val (@vals){
		    push @enumeration, {value=>$val}
		}
	    }
            if ( exists $selectable_labels{$label} || @enumeration) {
		if(@enumeration){
		    push @{$sections{$section_type}->[-1]->{parameters}}, { label => $label, 
                                                        pretty_label => $pretty_label, 
                                                        values => \@enumeration,
                                                        selectable => 1,
                                                        enumerated => 1,
                                                        section => $section,
                                                        comment => $component_cfg->get_comment_html($section, $parameter) };
		}
		else{
		    push @{$sections{$section_type}->[-1]->{parameters}}, { label => $label, 
									    pretty_label => $pretty_label, 
									    values => [{value=>$value}],
									    selectable => 1,
									    section => $section,
									    comment => $component_cfg->get_comment_html($section, $parameter) };            
		}
            } else {
                push @{$sections{$section_type}->[-1]->{parameters}}, { label => $label, 
                                                        pretty_label => $pretty_label, 
                                                        values => [{value=>$value}],
                                                        selectable => 0,
                                                        section => $section,
                                                        comment => $component_cfg->get_comment_html($section, $parameter) };
            }                
        }
    }
    
    $component_found++;
}

$tmpl->param( COMPONENT_FOUND => $component_found );
## $tmpl->param( COMPONENT_NAME    => $component_name );
$tmpl->param( BASIC_SECTIONS => $sections{basic} );
$tmpl->param( ADVANCED_SECTIONS => $sections{advanced} );
$tmpl->param( COMPONENT_ID => $component_id );

print $tmpl->output;

exit(0);
