#!/usr/bin/perl -w

use strict;
use CGI;
use CGI::Carp qw(fatalsToBrowser);
use Date::Manip;
use Ergatis::Common;
use Ergatis::ConfigFile;
use Ergatis::Monitor;
use HTML::Template;
use POSIX;
use XML::LibXML;

my $q = new CGI;

my $tmpl = HTML::Template->new( filename => 'templates/pipeline_list.tmpl',
                                die_on_bad_params => 0,
                                global_vars => 1
                              );


my $repository_root = $q->param("repository_root") || die "pass a repository root";
my $view_page_num = $q->param("page_num") || 1;
my $view = $q->param('view') || 'pipeline';
my $max_pages_to_show = 10;

## Quota information only availbe if project is in /usr/local/annotation/...
my $quotastring = &quota_string($repository_root);

my $ergatis_cfg = new Ergatis::ConfigFile( -file => "ergatis.ini" );
my $display_codebase = $ergatis_cfg->val( 'display_settings', 'display_codebase') || 0;

my $username = user_logged_in($ergatis_cfg);

my $pipeline_root = "$repository_root/workflow/runtime/pipeline";
my $project_conf_path = "$repository_root/workflow/project.config";
my $project_comment_file = "$repository_root/workflow/project.comment";
my $pipelines_per_page = $ergatis_cfg->val( 'display_settings', 'pipelines_per_page');
my $errors_found  = 0;
my $error_msgs = [];

## Should we limit our display list to just those pipelines that the user has created
my $per_account_pipelines = $ergatis_cfg->val('authentication', 'per_account_pipeline_security') || 0;

## Set the LOGGED_IN template variable if authentication is on and our user is logged in. This 
## variable controls whether or not the clone and archive/delete buttons are functional
my $auth_method = $ergatis_cfg->val('authentication', 'authentication_method');
my $logged_in = $auth_method eq 'open' || defined($username) ? 1: 0;

## Check user configured build access
my $require_login_for_pipeline_build = $ergatis_cfg->val('authentication', 'require_login_for_pipeline_build') || 0;

if (! -e $project_conf_path ) {
    &record_error("$project_conf_path not found");
    print $q->header( -type => 'text/html' );
    &print_template();
}

my $project_conf = new Ergatis::ConfigFile( -file => $project_conf_path );
my $ergatis_dir = $project_conf->val('project', '$;ERGATIS_DIR$;');

my $project_code = 0;
if ( $project_conf->val('project', '$;PROJECT_CODE$;') ) {
    $project_code = $project_conf->val('project', '$;PROJECT_CODE$;');
}

## Check for some required variables in the project.config
my @required_vars = ( '$;PROJECT$;', '$;REPOSITORY_ROOT$;', '$;TMP_DIR$;', '$;PROJECT_ID_REPOSITORY$;',
                      '$;ERGATIS_DIR$;', '$;LIB_DIR$;', '$;BIN_DIR$;', '$;DOCS_DIR$;' );
               
for my $var ( @required_vars ) {
    if (! defined $project_conf->val( 'project', $var ) ) {
        &record_error("$var not defined in project's configuration file");
    }
}

## Ensure that our required directories exist and are writable
my @required_dirs = ( $pipeline_root, "$repository_root/workflow", "$repository_root/workflow/runtime",
                      "$repository_root/workflow/lock_files", $project_conf->val( 'project', '$;PROJECT_ID_REPOSITORY$;' ) );

for my $dir ( @required_dirs ) {
    if (! -d $dir ) {
        &record_error("directory $dir not found");
    }
}

## If we have a project comment file grab its contents
my $project_comment = '';
if ( -e $project_comment_file ) {
    open(my $cfh, "<$project_comment_file") || die "can't read comment file: !";
    while (<$cfh>) {
        $project_comment .= $_;
    }
}

my $rdh;
if (! opendir $rdh, "$pipeline_root" ) {
    &record_error("can't open pipeline root directory: $!");
}

## If we've run into any errors at this point quit and print our template 
## with any errors that occurred.
if (scalar @{$error_msgs}) {
    print $q->header( -type => 'text/html' );
    &print_template();
}

## Do some groundwork to setup everything we need to paginate our pipelines.
my %pipeline_quickstats = get_pipeline_quickstats($ergatis_cfg, $rdh, $pipeline_root);
my $pipeline_count = scalar keys %pipeline_quickstats;

my @pagination_vars = prepare_pipeline_pagination($view_page_num, $repository_root, $max_pages_to_show, $pipelines_per_page, $pipeline_count);
my $page_count = $pagination_vars[0];
my $min_pipeline_pos = $pagination_vars[1];
my $max_pipeline_pos = $pagination_vars[2];
my $show_pre_continuation = $pagination_vars[3];
my $show_post_continuation = $pagination_vars[4];
my $page_links = $pagination_vars[5];

my @pipelines_to_parse = get_pipeline_ids_range($view, $min_pipeline_pos, $max_pipeline_pos, $pipelines_per_page, \%pipeline_quickstats);
my @pipelines_sorted = parse_pipeline_xmls($view, $pipeline_root, $repository_root, \%pipeline_quickstats, \@pipelines_to_parse);

## Populate the template
print $q->header( -type => 'text/html' );
print_template();


#############################################################
#                     SUBROUTINES                           #
#############################################################

###
# Render our template with corresponding variables
###
sub print_template {
    my $submenu_links = [ ];
                        
    if ( $require_login_for_pipeline_build && ! $username ) {
        push @$submenu_links, { label   => 'new pipeline (login)', 
                                is_last => 0, 
                                url     => "./login_form.cgi?redirect_url=" . CGI::escape("./build_pipeline.cgi?repository_root=$repository_root"),         
                              };
    } else {
        push @$submenu_links, { label   => 'new pipeline', 
                                is_last => 0, 
                                url     => "./build_pipeline.cgi?repository_root=$repository_root"
                              };
    }

    if ( $view eq 'component' ) {
        $tmpl->param( COMPONENTS => \@pipelines_sorted );
        push (@$submenu_links, { label   => 'view by pipeline', 
                                 is_last => 0, 
                                 url     => "./pipeline_list.cgi?repository_root=$repository_root" 
        });
        push (@$submenu_links, { label   => 'view by group', 
                                 is_last => 1, 
                                 url     => "./pipeline_list.cgi?repository_root=$repository_root&view=group" 
        });
    } elsif ( $view eq 'group' ) {
        $tmpl->param( PIPELINE_GROUPS => \@pipelines_sorted );
        push (@$submenu_links, { label   => 'view by pipeline', 
                                 is_last => 0, 
                                 url     => "./pipeline_list.cgi?repository_root=$repository_root" }
        );
        push (@$submenu_links, { label   => 'view by component', 
                                 is_last => 0, 
                                 url     => "./pipeline_list.cgi?repository_root=$repository_root&view=component" 
        });
        push (@$submenu_links, { label   => 'edit groups', 
                                 is_last => 1, 
                                 url     => "javascript:showGroupModificationMenu()" 
        });
    } else {
        $tmpl->param( PIPELINES => \@pipelines_sorted );
        push (@$submenu_links, { label   => 'view by component', 
                                 is_last => 0, 
                                 url     => "./pipeline_list.cgi?repository_root=$repository_root&view=component" 
        });
        push (@$submenu_links, { label   => 'view by group', 
                                 is_last => 1, 
                                 url     => "./pipeline_list.cgi?repository_root=$repository_root&view=group" 
        });
    }

    ## populate the template with the values that will always be passed.
    $tmpl->param( PROJECT_CODE     => $project_code );
    $tmpl->param( PROJECT_COMMENT  => $project_comment );
    $tmpl->param( PROJECT_COMMENT_FILE  => $project_comment_file );
    $tmpl->param( ERRORS_FOUND     => $errors_found );
    $tmpl->param( REPOSITORY_ROOT  => $repository_root );
    $tmpl->param( ERROR_MSGS       => $error_msgs );
    $tmpl->param( QUOTA_STRING     => $quotastring );
    $tmpl->param( ERGATIS_DIR      => $ergatis_dir );
    $tmpl->param( DISPLAY_CODEBASE => $display_codebase );
    $tmpl->param( COMPONENT_VIEW   => $view eq 'component' ? 1 : 0 );
    $tmpl->param( GROUP_VIEW       => $view eq 'group' ? 1 : 0 );
    $tmpl->param( QUICK_LINKS      => &get_quick_links($ergatis_cfg) );
    $tmpl->param( SUBMENU_LINKS    => $submenu_links );
    ## for pagination
    $tmpl->param( NEEDS_PAGINATION => $pipeline_count > $pipelines_per_page ? 1 : 0 );
    $tmpl->param( PIPELINE_COUNT   => $pipeline_count );
    if($page_links) {$tmpl->param( PAGE_LINKS       => $page_links )};
    $tmpl->param( IS_FIRST_PAGE    => $view_page_num == 1 ? 1 : 0 );
    $tmpl->param( IS_LAST_PAGE     => ($view_page_num * $pipelines_per_page) > $pipeline_count ? 1 : 0 );
    $tmpl->param( MIN_PIPELINE_POS => $min_pipeline_pos );
    $tmpl->param( MAX_PIPELINE_POS => $max_pipeline_pos );
    $tmpl->param( PREVIOUS_PAGE_URL => "./pipeline_list.cgi?repository_root=$repository_root&page_num=" . ($view_page_num - 1));
    $tmpl->param( NEXT_PAGE_URL    => "./pipeline_list.cgi?repository_root=$repository_root&page_num=" . ($view_page_num + 1));
    $tmpl->param( SHOW_PRE_CONTINUATION => $show_pre_continuation );
    $tmpl->param( SHOW_POST_CONTINUATION => $show_post_continuation );
    $tmpl->param( IS_PER_ACCOUNT_PIPELINES => $per_account_pipelines );
    $tmpl->param( LOGGED_IN => $logged_in);
    $tmpl->param( USER => $username );

    ## print the template
    print $tmpl->output();
    exit;
}


###
# Parses all pipeline ID's in the provided array and grab the corresponding 
# pipeline XML file that will be parsed and have its metadata stored in a 
# hash.
###
sub parse_pipeline_xmls {
    my ($view, $pipeline_root, $repository_root, $quickstats, $pipeline_ids_ref) = @_;
    my @pipeline_ids = @$pipeline_ids_ref;
    my $pipelines_metadata = {};


    my $parser = XML::LibXML->new();

    foreach my $pipeline_id (@pipeline_ids) {
        # Some variables to house metadata we want to track
        my $state = 'unknown';
        my $start_time = 'unknown';
        my $end_time = 'unknown';
        my $run_time = 'unknown';
        my $component_count = 0;
        my $component_aref = [];
        my $links_enabled = 1;
        my $error_message = 0;
        my $pipeline_comment = '';

        my $pipeline_file = "$pipeline_root/$pipeline_id/pipeline.xml"; 
        my $archive_link = "./archive_pipeline_form.cgi?repository_root=$repository_root&amp;pipeline_id=$pipeline_id";

        # Check if we have a lock file in place for this pipeline. 
        # A lockfile indicates that the pipeline could be in the archive/delete process
        my $lock_file = "$repository_root/workflow/lock_files/pipeline." .
                        "$pipeline_id.lock";
        if (-e $lock_file) {
            open (my $lockfh, "<$lock_file") || 
                die("Can't read lock file");
            $state = readline $lockfh;
            close ($lockfh);

            chomp $state;
        
            $run_time = ''; 
            $links_enabled = 0;
        } else {
            # If our pipeline file doesn't exists it might exist in gunzip'd
            # format.
            if (! -e $pipeline_file) {
                if (-e "$pipeline_file.gz") {
                    $pipeline_file = "$pipeline_file.gz";
                } else {
                    warn("Could not find pipeline XML for pipeline" .
                                  "$pipeline_id");
                    next;
                }
            }

            if (! -s $pipeline_file) {
                warn("Empty pipeline XML for pipeline $pipeline_id");
                next;
            }

            ## is there a comment?
            if ( -e "$pipeline_file.comment" ) {
                open(my $ifh, "<$pipeline_file.comment") || die "can't read comment file: $!";
                while (<$ifh>) {
                    $pipeline_comment .= $_;
                }
            }

            # Make sure that we actually have components to parse 
            my $doc = $parser->parse_file($pipeline_file);
            my $root = $doc->getDocumentElement;

            next if (! $doc->find('/commandSetRoot/commandSet/state') );
        
            my $commandSet = ($root->findnodes('commandSet'))[0];
            $state = $commandSet->findvalue('state') || "unknown";

            ## if the state is 'incomplete', check for a token file that indicates
            #   that this pipeline was submitted to a job manager.  this allows us to
            #   show a 'pending' state of the parent pipeline before the XML is parsed.
            if ( $state eq 'incomplete' && -e "$pipeline_file.submitted" ) {
                $state = 'pending';
            }

            ($start_time, $end_time, $run_time) = time_info_libxml($commandSet);

            # Build our component list
            if ($view eq 'component') {
                $component_aref = get_component_hash($commandSet);
            } else {
                ($component_count, $component_aref) = get_component_list($commandSet);
            }

            ## depending on the state, grab the top-level error message if there is one
            ##  there's no use doing this for running pipelines, since workflow won't
            ##  propagate the error until it's finished.
            if ( $state eq 'error' || $state eq 'failed' ) {
                $error_message = $commandSet->findvalue('status/message') || 0;
                $error_message = CGI::escapeHTML($error_message);
            }
            
            $pipelines_metadata->{$pipeline_id} = {
                'state'             => $state,
                'is_running'        => $state eq 'running' ? 1 : 0,
                'run_time'          => $run_time,
                'components'        => $component_aref,
                'component_count'   => $component_count,
                'component_label'   => $component_count > 1 ? ' components' : ' component',
                'view_link'         => "./view_pipeline.cgi?instance=$pipeline_file",
                'clone_link'        => "clone_pipeline.cgi?instance=$pipeline_file&repository_root=$repository_root",
                'archive_link'      => $archive_link,
                'links_enabled'     => $links_enabled,
                'error_message'     => $error_message,
                'has_comment'       => length($pipeline_comment),
                'pipeline_comment'  => $pipeline_comment,
                'pipeline_id'       => $pipeline_id,
                'last_mod'          => $quickstats->{$pipeline_id}->{'last_mod'},
                'pipeline_user'     => $quickstats->{$pipeline_id}->{'pipeline_user'}
            };
        }
    }

    return sort_pipelines($view, $pipelines_metadata);
}

###
# Sort our pipelines based off how the pipelines are being viewed:
#
#       pipeline - Show pipelines in chronological order
#       group - Chronological order by group
#       component - Sort by components
### 
sub sort_pipelines {
    my ($view, $pipelines_ref) = @_;
    my %pipelines = %$pipelines_ref;
    my @pipelines_sorted;
    my %pipeline_groups;
    my %components;

    if ( $view eq 'component' ) {
        ## Transformation into data structure passable to HTML::Template
        for my $pipeline ( %pipelines ) {
            next if (! exists $pipelines{$pipeline}{components});
        
            for my $component_ref ( @{$pipelines{$pipeline}{components}} ) {
                push (@{$components{ $$component_ref{name} }}, {
                    pipeline_id   => $pipelines{$pipeline}{pipeline_id},
                    pipeline_user => $pipelines{$pipeline}{pipeline_user},
                    pipeline_view_link => $pipelines{$pipeline}{view_link},
                    component_state => $$component_ref{state} || 'unknown',
                    component_token => $$component_ref{token} || 'unknown',
                    component_run_time => $$component_ref{run_time},
                    component_view_link => "./view_component.cgi?pipeline_xml=$repository_root" .
                                           "/workflow/runtime/$$component_ref{name}/" .
                                           "$pipelines{$pipeline}{pipeline_id}_" .
                                           "$$component_ref{token}/component.xml",
                });
                
                if ( -e "$repository_root/workflow/runtime/$$component_ref{name}/" .
                        "$pipelines{$pipeline}{pipeline_id}_$$component_ref{token}/" .
                        "$$component_ref{name}.$$component_ref{token}.final.config" ) {
                    $components{ $$component_ref{name} }[-1]{component_config_link} = 
                            "./view_formatted_ini_source.cgi?pipeline_id=$pipelines{$pipeline}{pipeline_id}&" .
                            "file=$repository_root/workflow/runtime/$$component_ref{name}/" .
                            "$pipelines{$pipeline}{pipeline_id}_$$component_ref{token}" .
                            "/$$component_ref{name}.$$component_ref{token}.final.config";
                } else {
                    $components{ $$component_ref{name} }[-1]{component_config_link} = 
                            "./view_formatted_ini_source.cgi?pipeline_id=$pipelines{$pipeline}{pipeline_id}&" .
                            "file=$repository_root/workflow/runtime/$$component_ref{name}/" .
                            "$pipelines{$pipeline}{pipeline_id}_$$component_ref{token}/" .
                            "$$component_ref{name}.$$component_ref{token}.user.config";
                }
            }
        }
        
        ## Second stage
        for ( sort keys %components ) {
            push (@pipelines_sorted, { name => $_, instances => $components{$_} });
        }

    } elsif ( $view eq 'group' ) {
        ## Put the pipelines in chronological order into groups
        for my $pipeline ( sort { $pipelines{$b}{last_mod} cmp $pipelines{$a}{last_mod} } keys %pipelines ) {
            $pipelines{$pipeline}{last_mod} = localtime( $pipelines{$pipeline}{last_mod} );
            
            my @groups = read_pipeline_groups( "$pipeline_root/$pipeline/pipeline.xml" );
            for my $group ( @groups ) {
                push (@{$pipeline_groups{$group}}, $pipelines{$pipeline});
                ${$pipeline_groups{$group}}[-1]{last_mod} = localtime( ${$pipeline_groups{$group}}[-1]{last_mod} );
            }
        }    

        for my $group ( keys %pipeline_groups ) {
            push (@pipelines_sorted, { name => $group, pipelines => \@{$pipeline_groups{$group}} });
        }

    } else {
        ## Just sort by straight chronological order
        for my $pipeline ( sort { $pipelines{$b}{last_mod} cmp $pipelines{$a}{last_mod} } keys %pipelines ) {
            push (@pipelines_sorted, $pipelines{$pipeline});
            $pipelines_sorted[-1]{last_mod} = localtime( $pipelines_sorted[-1]{last_mod} );
        }
    }

    return @pipelines_sorted;
}

###
# Reads pipeline groups given the path to a pipeline file.
###
sub read_pipeline_groups {
    my $pipeline_path = shift;
    my @groups = ();
    
    if ( -e "$pipeline_path.groups" ) {
        my $groupsfh = get_conditional_read_fh("$pipeline_path.groups");
        
        ## each line is a group
        while ( <$groupsfh>) {
            chomp;
            push @groups, $_;
        }
    }
    
    if ( scalar @groups < 1 ) {
        push @groups, 'ungrouped';
    }
    
    return @groups;
}

###
# Keeps track of an errors that may occur during the pipeline parsing
# process and stores them in an array to report to the user during 
# template-render execution
###
sub record_error {
    my $error_string = shift;
    
    $errors_found++;
    push @$error_msgs, { msg => $error_string };
}
