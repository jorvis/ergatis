#!/usr/local/bin/perl -w

use strict;
use CGI;
use CGI::Carp qw(fatalsToBrowser);
use Config::IniFiles;
use GO::Parser;
use HTML::Template;

my $q = new CGI;
print $q->header( -type => 'text/html' );

my $repository_root = $q->param('repository_root') || die "need a repository root";
my $shared_cfg = new Config::IniFiles( -file => "$repository_root/workflow_config_files/sharedconf.ini" );
my $workflowdocs_dir = $shared_cfg->val( 'init', '$;WORKFLOWDOCS_DIR$;' );

## path to the analysis ontology which contains the component list.
my $component_obo = 'analysis.obo';

my $tmpl = HTML::Template->new( filename => 'templates/creator.tmpl',
                                die_on_bad_params => 1,
                              );

################################
## build the component templates
## html definition lists would perhaps be most appropriate to hold the key/value pairs
## of the component configurations.
## 
my $component_template_html = '';
for my $conf_file ( glob "$workflowdocs_dir/*conf.ini" ) {
    $conf_file =~ m|.+\/(.+)conf.ini|;
    my $component_name = $1;
    
    $component_template_html .= "    <dl id='${component_name}_template'>\n";
    
    ## read this config file
    my $component_cfg = new Config::IniFiles( -file => $conf_file );
    
    ## it's possible that later the first dd could be comments and the second the value
    for my $section ( $component_cfg->Sections() ) {
        for my $parameter ( $component_cfg->Parameters($section) ) {
            $component_template_html .= "        <dt>$parameter</dt>\n";
            $component_template_html .= "        <dd>" . ($component_cfg->val($section, $parameter) || '') . "</dd>\n";
        }
    }
    
    $component_template_html .= "    </dl>\n";
}


##########################
## component chooser build
## we have to build this now and just push it out to the template since
##  HTML::Template couldn't handle the arbitrary levels of nesting.
my $component_chooser_html = '';

## parse the ontology file
my $parser = new GO::Parser( { handler => 'obj' } );
$parser->parse('analysis.obo');

## create a graph and pull the root where the component defs live
my $graph = $parser->handler->graph;
my $analysis_tool_root = $graph->get_term_by_name('analysis_tool');

## build the html from the ontology
&build_component_chooser_html( $analysis_tool_root, 3 );

## push the build out to the template
$tmpl->param( COMPONENT_CHOOSER_HTML  => $component_chooser_html );
$tmpl->param( COMPONENT_TEMPLATE_HTML => $component_template_html );
################

print $tmpl->output;

exit(0);

## recursive function used to build the component chooser menu
sub build_component_chooser_html {
    my ($parent, $indent_level) = @_;
    
    for my $node ( sort {$a->name cmp $b->name} @{ $graph->get_child_terms( $parent->acc ) } ) {
        ## perform any operations on the displayed name here
        my $name = $node->name();
        
        ## replace any underscores with spaces
        $name =~ s/_/ /g;
        
        ## if this node has any children, it is a sub-menu
        if ( $graph->n_children( $node->acc ) ) {
            add_to_component_chooser($indent_level, "<li class='menuparent'><a href='#'>$name</a>\n");
            add_to_component_chooser($indent_level, "    <ul>\n");
            
            build_component_chooser_html( $node, $indent_level + 2 );
            
            add_to_component_chooser($indent_level, "    </ul>\n");
            add_to_component_chooser($indent_level, "</li>\n");
            
        } else {
            add_to_component_chooser($indent_level, "<li><a href=\"javascript:addComponent('$name');\">$name</a></li>\n");
        }
    }
}

sub add_to_component_chooser {
    my ($l, $s) = @_;
    
    $component_chooser_html .= '    ' x $l;
    $component_chooser_html .= $s;
}
