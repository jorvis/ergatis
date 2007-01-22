#!/usr/local/bin/perl -w

use strict;
use CGI;
use CGI::Carp qw(fatalsToBrowser);
use Ergatis::Common;
use Ergatis::ConfigFile;
use HTML::Template;

my $q = new CGI;
print $q->header( -type => 'text/html' );

my $tmpl = HTML::Template->new( filename => 'templates/edit_config.tmpl',
                                die_on_bad_params => 0,
                                global_vars => 1
                              );

## read the ergatis config file
my $ergatis_cfg = new Ergatis::ConfigFile( -file => "ergatis.ini" );

my $config_path = $q->param('config') || die "config is a required parameter";
my $save = $q->param('save') || 0;
my $sections = [];

## make sure the configuration exists.
if ( -e $config_path ) {
    if (! -w $config_path ) {
        print_error_page( ergatis_cfg => $ergatis_cfg, 
                          message => "<p>The config file below is not writeable:</p><p>$config_path</p><p><a href='./edit_config.cgi?config=$config_path'>try again</a></p>",
                          links => [ ] 
                        );
        exit();
    }

    ## read this config file
    my $cfg = new Ergatis::ConfigFile( -file => $config_path );
    
    if ( $save ) {
    
        $cfg->import_form_data( $q );
        $cfg->RewriteConfig();
    
    } else {
    
        ## it's possible that later the first dd could be comments and the second the value
        for my $section ( $cfg->Sections() ) {

            push @$sections, { name => $section };

            for my $parameter ( $cfg->Parameters($section) ) {
                ## strip the parameter of the $; symbols
                my $label = $parameter;
                $label =~ s/\$\;//g;

                ## our variable naming 
                my $pretty_label = lc($label);
                   $pretty_label =~ s/_/ /g;

                push @{$$sections[-1]->{parameters}}, { label => $label, 
                                                        pretty_label => $pretty_label, 
                                                        value => $cfg->val($section, $parameter),
                                                        section => $section,
                                                        comment => $cfg->get_comment_html($section, $parameter) };
            }
        }
    }
} else {
    print_error_page( ergatis_cfg => $ergatis_cfg, 
                      message => "<p>The config file below does not exist:</p><p>$config_path</p><p><a href='./edit_config.cgi?config=$config_path'>try again</a></p>",
                      links => [ ] 
                    );
    exit();
}

$tmpl->param( SECTIONS            => $sections );
$tmpl->param( SAVE                => $save );
$tmpl->param( CONFIG              => $config_path );
$tmpl->param( QUICK_LINKS         => &get_quick_links($ergatis_cfg) );
$tmpl->param( SUBMENU_LINKS       => [
                                        #{ label => 'create project', is_last => 1, url => './create_project.cgi' },
                                     ] );
print $tmpl->output;

exit(0);
