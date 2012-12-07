#!/usr/bin/perl -w

=head1  SUMMARY

This script is used to edit any INI file within Ergatis.  Depending on the value of the
'save' parameter passed, it either presents a new form for user input or saves the form.

When save = 0 two additional parameters can be passed to make the script handle the 
special case of new project config files.  If you pass, for example:

    save = 0
    mode = new_project
    project_directory = /usr/local/projects/foo
    config=/usr/local/projects/foo/workflow/project.config

It will parse and display the config as a form but will replace values in many of the fields
based on default values set in the ergatis.ini

=cut

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

my $username = user_logged_in($ergatis_cfg);
my $auth_method = $ergatis_cfg->val('authentication', 'authentication_method');
unless ($auth_method eq 'open' || defined($username)) {
    print_error_page( ergatis_cfg => $ergatis_cfg,
                      message => "You must be logged in to edit project configs",
                      links => [],
                    );
    exit(0);
}

my $config_path = $q->param('config') || die "config is a required parameter";
my $save = $q->param('save') || 0;
my $sections = [];

## if mode = new_project the project_directory must also be passed
my $creating_new_project = 0;
if ( $q->param('mode') && $q->param('mode') eq 'new_project' ) {
    $creating_new_project = 1;
    
    if ( ! defined $q->param('project_directory') ) {
        print_error_page( ergatis_cfg => $ergatis_cfg, 
                          message => "<p>System error: attempted to create a new project config without " .
                                     "passing the project_directory parameter.  Please use the 'bugs' link " .
                                     "at the top of this page to submit a bug report</p>",
                          links => [ ] 
                        );
        exit();
    }
}

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
                
                ## if processing a new project form, we replace some of the default values here
                if ( $creating_new_project ) {
                    
                    ## the repository root value should be set to the passed project directory
                    if ( $section eq 'project' && $parameter eq '$;REPOSITORY_ROOT$;' ) {
                        $cfg->setval( $section, $parameter, $q->param('project_directory') );
                    
                    ## tmp dir should be the temp_space parameter from ergatis.ini with . '/$;PROJECT$;/$;COMPONENT_NAME$;/$;PIPELINEID$;_$;OUTPUT_TOKEN$;'
                    } elsif ( $section eq 'project' && $parameter eq '$;TMP_DIR$;' ) {
                        $cfg->setval( $section, $parameter, $ergatis_cfg->val('paths', 'temp_space') . '/$;PROJECT$;/$;COMPONENT_NAME$;/$;PIPELINEID$;_$;OUTPUT_TOKEN$;');
                    
                    ## the project_id_repository path should be standard from the project_dir
                    } elsif ( $section eq 'project' && $parameter eq '$;PROJECT_ID_REPOSITORY$;' ) {
                        $cfg->setval( $section, $parameter, $q->param('project_directory') . '/workflow/project_id_repository' );
                    
                    ## the next set should be based on the default_ergatis_dir path set in ergatis.ini
                    } elsif ( $section eq 'project' && $parameter eq '$;ERGATIS_DIR$;' ) {
                        $cfg->setval( $section, $parameter, $ergatis_cfg->val('paths', 'default_ergatis_dir')  );
                    
                    } elsif ( $section eq 'project' && $parameter eq '$;LIB_DIR$;' ) {
                        $cfg->setval( $section, $parameter, $ergatis_cfg->val('paths', 'default_ergatis_dir') . '/lib'  );
                    
                    } elsif ( $section eq 'project' && $parameter eq '$;BIN_DIR$;' ) {
                        $cfg->setval( $section, $parameter, $ergatis_cfg->val('paths', 'default_ergatis_dir') . '/bin'  );
                    
                    } elsif ( $section eq 'project' && $parameter eq '$;DOCS_DIR$;' ) {
                        $cfg->setval( $section, $parameter, $ergatis_cfg->val('paths', 'default_ergatis_dir') . '/docs'  );
                    
                    } elsif ( $section eq 'include' && $parameter eq '$;SOFTWARE_CONFIG$;' ) {
                        $cfg->setval( $section, $parameter, $ergatis_cfg->val('paths', 'default_ergatis_dir') . '/software.config'  );
                    }
                }
                
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
$tmpl->param( CREATING_NEW_PROJECT => $creating_new_project );
$tmpl->param( QUICK_LINKS         => &get_quick_links($ergatis_cfg) );
$tmpl->param( SUBMENU_LINKS       => [
                                        #{ label => 'create project', is_last => 1, url => './create_project.cgi' },
                                     ] );
print $tmpl->output;

exit(0);
