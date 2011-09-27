#!/usr/bin/perl -w

use strict;
use CGI;
use CGI::Carp qw(fatalsToBrowser);
use Ergatis::Common;
use Ergatis::ConfigFile;
use Ergatis::Monitor;
use HTML::Template;
use POSIX;
use XML::LibXML;

my $q = new CGI;

print $q->header( -type => 'text/html' );
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

my $pipeline_root = "$repository_root/workflow/runtime/pipeline";
my $project_conf_path = "$repository_root/workflow/project.config";
my $project_comment_file = "$repository_root/workflow/project.comment";
my $pipelines_per_page = $ergatis_cfg->val( 'display_settings', 'pipelines_per_page');
my $errors_found  = 0;
my $error_msgs = [];

if (! -e $project_conf_path ) {
    &record_error("$project_conf_path not found");
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
    } else {
        ## make sure it is writeable
        if (! -w $dir) {
            &record_error("$dir not writable");
        }
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
&print_template() if scalar @{$error_msgs};

## Do some groundwork to setup everything we need to paginate our pipelines.
my %pipeline_quickstats = get_pipeline_quickstats($rdh);
my $pipeline_count = scalar keys %pipeline_quickstats;

my @pagination_vars = prepare_pipeline_pagination($view_page_num, $max_pages_to_show, $pipelines_per_page, $pipeline_count);
my $page_count = $pagination_vars[0];
my $min_pipeline_pos = $pagination_vars[1];
my $max_pipeline_pos = $pagination_vars[2];
my $show_pre_continuation = $pagination_vars[3];
my $show_post_continuation = $pagination_vars[4];
my $page_links = $pagination_vars[5];

my @pipelines_to_parse = get_pipeline_ids_range($view, $min_pipeline_pos, $max_pipeline_pos, \%pipeline_quickstats);
my @pipelines_sorted = parse_pipeline_xmls($view, $pipeline_root, $repository_root, \%pipeline_quickstats, \@pipelines_to_parse);

## Populate the template
print_template();

#############################################################
#                     SUBROUTINES                           #
#############################################################

###
# Render our template with corresponding variables
###
sub print_template {
    my $submenu_links = [
                            { label   => 'new pipeline', 
                              is_last => 0, 
                              url     => "./build_pipeline.cgi?repository_root=$repository_root"
                            },
                        ];

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
# Generates a hash of all the components 
###
sub get_component_hash {
    my $command = shift;
    my @components;

    foreach my $node ($command->findnodes('//commandSet')) {
        my $name = $node->findvalue('name');
        if ($name && $name =~ /^(.+?)\.(.+)/) {
            push (@components, { name => $1, token => $2 });

            my $state = $node->findvalue('state');
            if ($state) {
                $components[-1]{'state'} = $state;

                if ($state eq 'error' || $state eq 'failed') {
                    $components[-1]{error_count}++;
                }
            }
            
            ($components[-1]{'start_time'}, $components[-1]{'end_time'}, $components[-1]{'run_time'}) = time_info_libxml($node);
        }
    }

    return \@components;
}

###
# Generates a list of all components in this pipeline and a count of each 
# component.
###
sub get_component_list {
    my $command = shift;
    my %components;
    my $component_aref = [];
    my $component_count = 0;

    foreach my $node ($command->findnodes('//commandSet')) {
        my $name = $node->findvalue('name');
        if ($name && $name=~ /^(.+?)\./) {
            $components{$1}{'count'}++;
    
            my $state = $node->findvalue('state');
            if ($state) {
                if ($state eq 'error' || $state eq 'failed') {
                    $components{$1}{'error_count'}++;
                }
            }
        }
    }

    # HTML::Template requires us to pass any complex data structures in an 
    # array of hashes so we need to convert
    foreach my $component (sort keys %components) {
        $component_count++;
        push (@{ $component_aref }, {
            'name'          => $component,
            'count'         => $components{$component}{'count'},
            'error_count'   => $components{$component}{'error_count'}
        });
    }

    return ($component_count, $component_aref);
}

###
# Grabs the relevant time information from our pipeline XML
###
sub time_info_libxml {
    my $command = shift;

    my $start_time_str = $command->findvalue('startTime');
    if (! $start_time_str ) {
        return ('unavailable', 'unavailable', 'unavailable');
    }
    
    my $start_time_obj = ParseDate($start_time_str);
    my $start_time = UnixDate($start_time_obj, "%c");

    my ($end_time_obj, $end_time);
    my $end_time_str = $command->findvalue('endTime');
    if ($end_time_str) {
        $end_time_obj = ParseDate($end_time_str);
        $end_time = UnixDate($end_time_obj, "%c");
    }

    my $runtime = "?";

    if ($start_time_obj) {
        my $diffstring;
        $runtime = '';
        
        if ($end_time_obj) {
            $diffstring = DateCalc($start_time_obj, $end_time_obj);
        } else {
            $diffstring = DateCalc($start_time_obj, "now");
        }
        
       ## take out any non \d: characters
        $diffstring =~ s/[^0-9\:]//g;

        my @parts = split(/:/, $diffstring);
        
        ## years + months + weeks + days
        my $val = ($parts[0] * 365) + ($parts[1] * 30) + ($parts[2] * 7) + ($parts[3]);
        if ($val > 1) {
            $runtime .= "$val days ";
        } elsif ($val == 1) {
            $runtime .= "$val day ";
        }
        
        $runtime .= "$parts[4] hr " if $parts[4];
        $runtime .= "$parts[5] min " if $parts[5];
        $runtime .= "$parts[6] sec";
    }

    $runtime = '&lt; 1 sec' if $runtime eq '0 sec';

    return ($start_time, $end_time, $runtime);
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
                            "./view_formatted_ini_source.cgi?file=$repository_root" .
                            "/workflow/runtime/$$component_ref{name}/" .
                            "$pipelines{$pipeline}{pipeline_id}_$$component_ref{token}" .
                            "/$$component_ref{name}.$$component_ref{token}.final.config";
                } else {
                    $components{ $$component_ref{name} }[-1]{component_config_link} = 
                            "./view_formatted_ini_source.cgi?file=$repository_root/" .
                            "workflow/runtime/$$component_ref{name}/" .
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
# Generates the range of pipeline ID's that should be displayed for this
# specific page of pipelines.
###
sub get_pipeline_ids_range {
    my ($view, $min_pipeline_pos, $max_pipeline_pos, $quickstats_ref) = @_;
    my %pipeline_quickstats = %$quickstats_ref;
    my @pipelines_to_parse = ();

    if ( $view eq 'group' || $view eq 'component' ) {
        ## No current pagination on grouped or component-based pipeline list views
        @pipelines_to_parse = keys %pipeline_quickstats;
    } else {
        ## Sort and record which to parse
        my $current_pipeline_pos = 0;
        for my $pipeline_id ( sort { $pipeline_quickstats{$b}{last_mod} cmp $pipeline_quickstats{$a}{last_mod} } keys %pipeline_quickstats ) {
            $current_pipeline_pos++;
            
            ## Check and see if this pipeline is in our page range
            next if ( $current_pipeline_pos < $min_pipeline_pos || $current_pipeline_pos > $max_pipeline_pos );
            
            push (@pipelines_to_parse, $pipeline_id);
            
            if ( scalar(@pipelines_to_parse) >= $pipelines_per_page ) {
                last;
            }
        }
    }

    return @pipelines_to_parse;
}

###
# Calculates several statistics/values that are needed to properly pagination 
# the pipelines for this project.
###
sub prepare_pipeline_pagination {
    my ($view_page_num, $max_pages_to_show, $pipelines_per_page, $pipeline_count) = @_;
    my $page_count;
    
    if ( $pipeline_count > $pipelines_per_page ) {
        if ( $pipeline_count % $pipelines_per_page ) {
            $page_count = int($pipeline_count / $pipelines_per_page) + 1;
        } else {
            $page_count = int($pipeline_count / $pipelines_per_page);
        }
    } else {
        $page_count = 1;
    }

    ## Calculate the min/max pipeline positions to show based on the page requested
    my $min_pipeline_pos;
    my $max_pipeline_pos;

    if ( $view_page_num == 1 ) {
        $min_pipeline_pos = 1;
    } else {
        $min_pipeline_pos = ($pipelines_per_page * ($view_page_num - 1)) + 1;
    }

    if ( $view_page_num * $pipelines_per_page >= $pipeline_count ) {
        $max_pipeline_pos = $pipeline_count;
    } else {
        $max_pipeline_pos = $view_page_num * $pipelines_per_page;
    }

    my $page_links = [];

    ## We can't show ALL pages since there could be dozens.  
    ## Calculate the next highest/lowest interval of '10', then hide anything outside of that range
    my ($next_lower_interval, $next_higher_interval);

    if ( $view_page_num == $max_pages_to_show ) {
        $next_lower_interval  = ( $view_page_num / $max_pages_to_show );
    } elsif ( $view_page_num % $max_pages_to_show ) {
        $next_lower_interval  = floor($view_page_num/$max_pages_to_show)*$max_pages_to_show + 1;
    } else {
        $next_lower_interval  = $view_page_num - $max_pages_to_show + 1;
    }

    $next_higher_interval = ceil($view_page_num/$max_pages_to_show)*$max_pages_to_show + 1;

    ## These are calculated later and control the display of the boxes with '...' on either
    ## of the pagination ranges
    my $show_pre_continuation  = $next_lower_interval > 1 ? 1 : 0;
    my $show_post_continuation = $next_higher_interval < $page_count ? 1 : 0;

    for ( my $i=1; $i<=$page_count; $i++ ) {
        next if $i >= $next_higher_interval || $i < $next_lower_interval;

        push @$page_links, {
            page_num => $i,
            is_active => $i == $view_page_num ? 1 : 0,
            url => "./pipeline_list.cgi?repository_root=$repository_root&page_num=$i"
        };
    }

    return  ($page_count, $min_pipeline_pos, $max_pipeline_pos, $show_pre_continuation, 
             $show_post_continuation, $page_links);
}

###
# Retrieves all pipeline XML for the given project and some preliminary metadata 
# on each pipeline to aid in pagination.
###
sub get_pipeline_quickstats {
    my $rdh = shift;
    my %pipeline_quickstats;

    foreach my $pipeline_id ( readdir $rdh ) {
        next unless ( $pipeline_id =~ /^\d+$/ );

        my $pipeline_file = "$pipeline_root/$pipeline_id/pipeline.xml";  ## may be modified below
        
        ## if only the pipeline.xml exists, we can do less
        if (! -e $pipeline_file ) {

            if ( -e "$pipeline_file.gz" ) {
                $pipeline_file .= '.gz';
            } else {
                ## didn't actually find a pipeline file
                next;
            }
        }
       
        my $filestat = stat($pipeline_file);
        my $pipeline_user = getpwuid($filestat->uid);
        my $last_mod = $filestat->mtime;
       
        $pipeline_quickstats{$pipeline_id} = { 
                            pipeline_id     => $pipeline_id,  ## set again so we can directly export to HTML::Template
                            last_mod        => $last_mod,
                            path            => $pipeline_file,
                            pipeline_user   => $pipeline_user,
        };
    }

    return %pipeline_quickstats;
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
