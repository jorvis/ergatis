#!/usr/local/bin/perl -w

use strict;
use CGI;
use CGI::Carp qw(fatalsToBrowser);
use Ergatis::Common;
use Ergatis::ConfigFile;
use HTML::Template;

my $q = new CGI;
print $q->header( -type => 'text/html' );

my $tmpl = HTML::Template->new( filename => 'templates/documentation.tmpl',
                                die_on_bad_params => 0,
                              );

## read the ergatis config file
my $ergatis_cfg = new Ergatis::ConfigFile( -file => 'ergatis.ini' );

my $article = $q->param('article') || 'main';
my $page = $q->param('page') || 'index';
my $ergatis_dir = $ergatis_cfg->val('paths', 'default_ergatis_dir' );
my $doc_path;

## either getting a component documentation list, or 
if ( $article eq 'components' && $page ne 'index' ) {
    $doc_path = "$ergatis_dir/docs/documentation/$page.tmpl";
    
} else {
    $doc_path = "templates/documentation/$article/$page.tmpl";
}

my $page_tmpl;

if ( -e $doc_path ) {
    $page_tmpl = HTML::Template->new( filename => $doc_path,
                                      die_on_bad_params => 0 );
                                      
    if ( $article eq 'components' && $page eq 'index' ) {
        my $component_sections = generate_component_list();
        $page_tmpl->param( COMPONENT_SECTIONS => $component_sections );
        $page_tmpl->param( SECTION_COUNT => scalar @$component_sections );
        $page_tmpl->param( ERGATIS_DIR => $ergatis_dir );
    }
    
} else {
    $page_tmpl = HTML::Template->new( filename => 'templates/documentation/main/not_found.tmpl',
                                      die_on_bad_params => 0 );
    $page_tmpl->param( ARTICLE => $article );
    $page_tmpl->param( ARTICLE_DISPLAY => display_version($article) );
    $page_tmpl->param( PAGE => $page );
    $page_tmpl->param( PAGE_DISPLAY => display_version($page) );
    $page_tmpl->param( DOC_PATH => $doc_path );
}


## populate the parent template
$tmpl->param( ARTICLE => $article );
$tmpl->param( ARTICLE_DISPLAY => display_version($article) );
$tmpl->param( PAGE => $page );
$tmpl->param( PAGE_DISPLAY => display_version($page) );
$tmpl->param( IS_INDEX => $page eq 'index' ? 1 : 0 );
$tmpl->param( DOC_CONTENT => $page_tmpl->output() );

$tmpl->param( QUICK_LINKS         => &get_quick_links($ergatis_cfg) );
$tmpl->param( SUBMENU_LINKS       => [
                                        #{ label => 'create project', is_last => 1, url => './create_project_form.cgi' },
                                     ] );

print $tmpl->output;


sub display_version {
    my $str = shift;
    $str =~ s/_/ /g;
    return $str;
}

sub generate_component_list {
    my %sections;

    opendir(my $idh, "$ergatis_dir/docs") || die "can't read component directory ($ergatis_dir/docs): $!";
    while ( my $thing = readdir $idh ) {
        if ( $thing =~ /(.+).config$/ ) {
            my $component_name = $1;
            my $status = $ergatis_cfg->component_status( $component_name );
            
            if (! -e "$ergatis_dir/docs/$thing" ) {
                die("$ergatis_dir/docs/$thing");
            }
            
            my $component_cfg = new Ergatis::ConfigFile( -file => "$ergatis_dir/docs/$thing" );
            my $class_string = 'unclassified';
            
            if ( $component_cfg && $component_cfg->val( 'interface', 'classification' ) ) {
                $class_string = $component_cfg->val( 'interface', 'classification' );
            }
            
            ## this can be a comma-separated list.  split it up
            for my $class ( split(',', $class_string) ) {
                ## strip leading/tailing whitespace
                $class =~ s/^\s+//g;
                $class =~ s/\s+$//g;
                
                my $documented = 0;
                
                if ( -e "$ergatis_dir/docs/documentation/$component_name.tmpl" ) {
                    $documented = 1;
                }
                    
                push @{$sections{$class}}, { name => $component_name, 
                                             documented => $documented, 
                                             disabled => $status eq 'disabled' ? 1 : 0
                                           };
            }
        }
    }

    my $return_sections = [];

    foreach my $section ( sort keys %sections ) {
        push @$return_sections, { section => $section, components => $sections{$section} };
    }

    return $return_sections;
}






