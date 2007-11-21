#!/usr/bin/perl -w

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
        my $doc_class = 'categorized';
        if ( $q->param('categorized') && $q->param('categorized') == 0 ) {
            $doc_class = 'alphabetical';
        }
    
        my $component_sections = generate_component_list( $doc_class );
        $page_tmpl->param( COMPONENT_GROUPS => $component_sections );
        $page_tmpl->param( SECTION_COUNT => scalar @$component_sections );
        $page_tmpl->param( ERGATIS_DIR => $ergatis_dir );
        $page_tmpl->param( DOC_CLASS => 'wha?');
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
    ## usually 'alphabetical' or 'categorized'
    my $display_method = shift;
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
            my $class_string = '';
            
            if ( $display_method eq 'categorized' ) {
                if ( $component_cfg && $component_cfg->val( 'interface', 'classification' ) ) {
                    $class_string = $component_cfg->val( 'interface', 'classification' );
                } else {
                    $class_string = 'unclassified';
                }
            } else {
                if ( $component_name =~ /^[A-M]/i ) {
                    $class_string = 'A-M';
                } elsif ( $component_name =~ /[N-Z]/i ) {
                    $class_string = 'N-Z';
                }
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
    ## 3 columns (groups)
    ## if you make changes to the group count you'll need to reflect them in the documentation.css
    my $component_groups = [ { group_num => 0, component_sections => [], display_method => $display_method },
                             { group_num => 1, component_sections => [], display_method => $display_method },
                             { group_num => 2, component_sections => [], display_method => $display_method },
                             { group_num => 3, component_sections => [] },
                           ];
    my $group_counts = [0,0,0,0];
    
    foreach my $section ( sort keys %sections ) {
        my $min_group_num = min_group($group_counts);
        
        ## remember how many are in this column (add a few to make up for the header overhead)
        $$group_counts[$min_group_num] += scalar @{$sections{$section}} + 2;
        
        ## sort the components within each section
        my @sorted_components = sort {$$a{name} cmp $$b{name}} @{$sections{$section}};
        
        push @{$$component_groups[$min_group_num]{component_sections}}, { section => $section, components => \@sorted_components };
    }

    return $component_groups;
}


sub min_group {
    my $counts = shift;
    my $min_num = 0;
    my $min_count = 10000000000000;
    
    for ( my $i=0; $i< scalar @$counts; $i++ ) {
        
        if ( $$counts[$i] < $min_count ) {
            $min_count = $$counts[$i];
            $min_num = $i;
        }
    }
    
    return $min_num;
}














