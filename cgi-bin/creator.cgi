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
my $shared_cfg = new Ergatis::ConfigFile( -file => "$repository_root/workflow/project.config" );
my $workflowdocs_dir = $shared_cfg->val( 'init', '$;WORKFLOWDOCS_DIR$;' );

## path to the analysis ontology which contains the component list.
my $component_obo = 'analysis.obo';

my $tmpl = HTML::Template->new( filename => 'templates/creator.tmpl',
                                die_on_bad_params => 1,
                              );


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
$tmpl->param( REPOSITORY_ROOT         => $repository_root );
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
