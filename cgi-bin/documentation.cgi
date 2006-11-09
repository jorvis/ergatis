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
my $ergatis_cfg = new Ergatis::ConfigFile( -file => "ergatis.ini" );

my $article = $q->param('article') || 'main';
my $page = $q->param('page') || 'index';
my $ergatis_dir = $ergatis_cfg->val('paths', 'default_ergatis_dir' );
my $doc_path;

if ( $article eq 'components' ) {
    $doc_path = "$ergatis_dir/docs/documentation/$page.tmpl";
} else {
    $doc_path = "templates/documentation/$article/$page.tmpl";
}

my $page_tmpl;

if ( -e $doc_path ) {
    $page_tmpl = HTML::Template->new( filename => $doc_path,
                                      die_on_bad_params => 0 );
} else {
    $page_tmpl = HTML::Template->new( filename => 'templates/documentation/main/not_found.tmpl',
                                      die_on_bad_params => 0 );
    $page_tmpl->param( ARTICLE => $article );
    $page_tmpl->param( ARTICLE_DISPLAY => display_version($article) );
    $page_tmpl->param( PAGE => $page );
    $page_tmpl->param( PAGE_DISPLAY => display_version($page) );
    $page_tmpl->param( DOC_PATH = $doc_path );
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
