#!/usr/bin/perl -w

use strict;
use CGI;
use CGI::Carp qw(fatalsToBrowser);
use Ergatis::Common;
use Ergatis::ConfigFile;
use Ergatis::Monitor;
use HTML::Template;
use POSIX;
use XML::Twig;

my $q = new CGI;

print $q->header( -type => 'text/html' );
my $tmpl = HTML::Template->new( filename => 'templates/pipeline_list.tmpl',
                                die_on_bad_params => 0,
                                global_vars => 1
                              );

my $repository_root = $q->param("repository_root") || die "pass a repository root";
my $view_page_num = $q->param("page_num") || 1;
my $max_pages_to_show = 10;

## quota information (only works if in /usr/local/annotation/SOMETHING)
my $quotastring = &quota_string($repository_root);

my $ergatis_cfg = new Ergatis::ConfigFile( -file => "ergatis.ini" );
my $display_codebase = $ergatis_cfg->val( 'display_settings', 'display_codebase') || 0;

my $pipeline_root = "$repository_root/workflow/runtime/pipeline";
my $project_conf_path = "$repository_root/workflow/project.config";
my $pipelines_per_page = $ergatis_cfg->val( 'display_settings', 'pipelines_per_page');
my $errors_found  = 0;
my $error_msgs = [];
my %components;

## make sure it has a project conf file, else throw the warning and quit
if (! -e $project_conf_path ) {
    &record_error("$project_conf_path not found");
    &print_template();
}

## pull the ergatis dir from the shared conf file
my $project_conf = new Ergatis::ConfigFile( -file => $project_conf_path );

## see if a project code is defined
my $project_code = 0;
if ( $project_conf->val('project', '$;PROJECT_CODE$;') ) {
    $project_code = $project_conf->val('project', '$;PROJECT_CODE$;');
}

#################
## check for some required variables in the project.config
my @required_vars = ( '$;PROJECT$;', '$;REPOSITORY_ROOT$;', '$;TMP_DIR$;', '$;PROJECT_ID_REPOSITORY$;',
                      '$;ERGATIS_DIR$;', '$;LIB_DIR$;', '$;BIN_DIR$;', '$;DOCS_DIR$;' );
               
for my $var ( @required_vars ) {
    if (! defined $project_conf->val( 'project', $var ) ) {
        &record_error("$var not defined in project's configuration file");
    }
}

#################
## a check a selection of required, writeable directories
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

my $rdh;
if (! opendir $rdh, "$pipeline_root" ) {
    &record_error("can't open pipeline root directory: $!");
}

## quit and throw the template if there were errors.
&print_template() if scalar @{$error_msgs};

my $ergatis_dir = $project_conf->val('project', '$;ERGATIS_DIR$;');

## views are either component or pipeline (how the list is organized)
my $view = $q->param('view') || 'pipeline';

## read through the directory just to get a quick list of what pipelines
#   are there along with their ages.  this enables faster paginated views
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

my $pipeline_count = scalar keys %pipeline_quickstats;
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

## calculate the min/max pipeline positions to show based on the page requested
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

## we can't show ALL pages since there could be dozens.  
## calculate the next highest/lowest interval of '10', then hide anything outside of that range
my ($next_lower_interval, $next_higher_interval);

if ( $view_page_num == $max_pages_to_show ) {
    $next_lower_interval  = ( $view_page_num / $max_pages_to_show );
} elsif ( $view_page_num % $max_pages_to_show ) {
    $next_lower_interval  = floor($view_page_num/$max_pages_to_show)*$max_pages_to_show + 1;
} else {
    $next_lower_interval  = $view_page_num - $max_pages_to_show + 1;
}

$next_higher_interval = ceil($view_page_num/$max_pages_to_show)*$max_pages_to_show + 1;

## these are calculated later and control the display of the boxes with '...' on either
#   of the pagination ranges
my $show_pre_continuation  = $next_lower_interval > 1 ? 1 : 0;
my $show_post_continuation = $next_higher_interval < $page_count ? 1 : 0;

for ( my $i=1; $i<=$page_count; $i++ ) {
    next if $i >= $next_higher_interval || $i < $next_lower_interval;

    push @$page_links, {
        page_num => $i,
        is_active => $i == $view_page_num ? 1 : 0,
        url => "./pipeline_list.cgi?repository_root=$repository_root&page_num=$i"
    };
    
    #last if scalar(@$page_links) == $max_pages_to_show;
}

## filter the full %pipelines has to include only those that need to be parsed.
#   create @ids to keep
#   populate a new %pipelines as we parse, TODO: rename old %pipelines
my @pipelines_to_parse = ();
my %pipelines;

if ( $view eq 'group' || $view eq 'component' ) {
    ## no current pagination on grouped or component-based pipeline list views
    @pipelines_to_parse = keys %pipeline_quickstats;
} else {
    ## sort and record which to parse
    my $current_pipeline_pos = 0;
    for my $pipeline_id ( sort { $pipeline_quickstats{$b}{last_mod} cmp $pipeline_quickstats{$a}{last_mod} } keys %pipeline_quickstats ) {
        $current_pipeline_pos++;
        
        ## check and see if this pipeline is in our page range
        next if ( $current_pipeline_pos < $min_pipeline_pos || $current_pipeline_pos > $max_pipeline_pos );
        
        push @pipelines_to_parse, $pipeline_id;
        
        if ( scalar(@pipelines_to_parse) >= $pipelines_per_page ) {
            last;
        }
    }
}

foreach my $pipeline_id ( @pipelines_to_parse ) {

    my $state = 'unknown';
    my $start_time = 'unknown';
    my $end_time = 'unknown';
    my $run_time = 'unknown';
    my $component_count = 0;
    my $component_aref = [];
    my $component_label = ' components';
    my $pipeline_file = "$pipeline_root/$pipeline_id/pipeline.xml";  ## may be modified below
    my $archive_link = "./archive_pipeline_form.cgi?repository_root=$repository_root&amp;pipeline_id=$pipeline_id";
    my $links_enabled = 1;
    my $error_message = 0;
    my $pipeline_comment = '';
    
    ## see if this pipeline is locked (usually for maintenance such as deleting or archiving.)
    my $lock_file = "$repository_root/workflow/lock_files/pipeline.$pipeline_id.lock";
    if ( -e $lock_file ) {
        ## get the state from the lock file, currently the only contents
        open (my $lockfh, "<$lock_file") || die "can't read lock file\n";
        $state = readline $lockfh;
        chomp $state;
    
        $run_time = '';
        $links_enabled = 0;
    
    } else {

        ## is there a comment?
        if ( -e "$pipeline_file.comment" ) {
            open(my $ifh, "<$pipeline_file.comment") || die "can't read comment file: $!";
            while (<$ifh>) {
                $pipeline_comment .= $_;
            }
        }

        ## if only the pipeline.xml exists, we can do less
        if (! -e $pipeline_file ) {

            if ( -e "$pipeline_file.gz" ) {
                $pipeline_file .= '.gz';
            } else {
                ## didn't actually find a pipeline file
                $pipeline_count--;
                next;
            }
        }

        if (! -s $pipeline_file ) {
            print STDERR "skipped empty pipeline file $pipeline_file\n";
            next;
        }

        my $ifh;
        if ($pipeline_file =~ /\.gz/) {
            open($ifh, "<:gzip", "$pipeline_file") || die "can't read $pipeline_file: $!"; 
        } else {
            open($ifh, "<$pipeline_file") || die "can't read $pipeline_file: $!";       
        }

        my $twig = new XML::Twig;
        $twig->parse($ifh);

        my $commandSetRoot = $twig->root;
        my $commandSet = $commandSetRoot->first_child('commandSet');

        next if (! $commandSet );

        if ( $commandSet->has_child('state') ) {
            $state  = $commandSet->first_child('state')->text();
        }

        ($start_time, $end_time, $run_time) = &time_info( $commandSet );
        
        ## depending on the state, grab the top-level error message if there is one
        ##  there's no use doing this for running pipelines, since workflow won't
        ##  propagate the error until it's finished.
        if ( $state eq 'error' || $state eq 'failed' ) {
            if ( $commandSet->has_child('status') && 
                 $commandSet->first_child('status')->has_child('message') ) {
            
                $error_message = $commandSet->first_child('status')->first_child('message')->text();

                ## handle illegal characters
                $error_message = CGI::escapeHTML($error_message);
            }
        }
        
        ## this is done as a new twig parse since elements can be nested
        ## at any level.
        if ( $view eq 'component' ) {
            $component_aref = &component_info_aref( $pipeline_file );

        } else {
        
            my %pipeline_components = &component_count_hash( $pipeline_file );

            foreach my $component (sort keys %pipeline_components) {
                $component_count += $pipeline_components{$component}{count};
                push @$component_aref, { name => $component, 
                                         count => $pipeline_components{$component}{count},
                                         error_count => $pipeline_components{$component}{error_count},
                                       };
            }

            ## reformat component_count to include a label
            if ($component_count == 1) {
                $component_label = ' component';
            }
        }
    }
    
    my $view_link = "./view_pipeline.cgi?instance=$pipeline_file";
    
    $pipelines{$pipeline_id} = { 
                        state           => $state,
                        is_running      => $state eq 'running' ? 1 : 0,
                        run_time        => $run_time,
                        components      => $component_aref,
                        component_count => $component_count,
                        component_label => $component_label,
                        view_link       => $view_link,
                        archive_link    => $archive_link,
                        clone_link      => "clone_pipeline.cgi?instance=$pipeline_file&repository_root=$repository_root",
                        links_enabled   => $links_enabled,
                        error_message   => $error_message,
                        has_comment     => length($pipeline_comment),
                        pipeline_comment => $pipeline_comment,
                        
                        ## these were set previously, but carried over here
                        pipeline_id     => $pipeline_id,
                        last_mod        => $pipeline_quickstats{$pipeline_id}{last_mod},
                        pipeline_user   => $pipeline_quickstats{$pipeline_id}{pipeline_user},
                      };
}

## sort the pipelines
my @pipelines_sorted;
my @components_sorted;
my %pipeline_groups;
my @pipeline_groups_sorted;

if ( $view eq 'component' ) {
    
    ## transformation into data structure passable to HTML::Template
    for my $pipeline ( %pipelines ) {
        next if (! exists $pipelines{$pipeline}{components});
    
        for my $component_ref ( @{$pipelines{$pipeline}{components}} ) {

            push @{$components{ $$component_ref{name} }}, {
                pipeline_id   => $pipelines{$pipeline}{pipeline_id},
                pipeline_user => $pipelines{$pipeline}{pipeline_user},
                pipeline_view_link => $pipelines{$pipeline}{view_link},
                component_state => $$component_ref{state} || 'unknown',
                component_token => $$component_ref{token} || 'unknown',
                component_run_time => $$component_ref{run_time},
                 component_view_link => "./view_component.cgi?pipeline_xml=$repository_root/workflow/runtime/$$component_ref{name}/$pipelines{$pipeline}{pipeline_id}_$$component_ref{token}/component.xml",
            };
            
            if ( -e "$repository_root/workflow/runtime/$$component_ref{name}/$pipelines{$pipeline}{pipeline_id}_$$component_ref{token}/$$component_ref{name}.$$component_ref{token}.final.config" ) {
                $components{ $$component_ref{name} }[-1]{component_config_link} = "./view_formatted_ini_source.cgi?file=$repository_root/workflow/runtime/$$component_ref{name}/$pipelines{$pipeline}{pipeline_id}_$$component_ref{token}/$$component_ref{name}.$$component_ref{token}.final.config";
            } else {
                $components{ $$component_ref{name} }[-1]{component_config_link} = "./view_formatted_ini_source.cgi?file=$repository_root/workflow/runtime/$$component_ref{name}/$pipelines{$pipeline}{pipeline_id}_$$component_ref{token}/$$component_ref{name}.$$component_ref{token}.user.config";
            }
        }
    }
    
    ## second stage
    for ( sort keys %components ) {
        push @components_sorted, { name => $_, instances => $components{$_} };
    }

} elsif ( $view eq 'group' ) {

    ## put the pipelines in chronological order into groups
    for my $pipeline ( sort { $pipelines{$b}{last_mod} cmp $pipelines{$a}{last_mod} } keys %pipelines ) {
        $pipelines{$pipeline}{last_mod} = localtime( $pipelines{$pipeline}{last_mod} );
        
        my @groups = read_pipeline_groups( "$pipeline_root/$pipeline/pipeline.xml" );
        
        for my $group ( @groups ) {
            push @{$pipeline_groups{$group}}, $pipelines{$pipeline};
            ${$pipeline_groups{$group}}[-1]{last_mod} = localtime( ${$pipeline_groups{$group}}[-1]{last_mod} );
        }
    }    

    
    for my $group ( keys %pipeline_groups ) {
        push @pipeline_groups_sorted, { name => $group, pipelines => \@{$pipeline_groups{$group}} };
    }

} else {

    for my $pipeline ( sort { $pipelines{$b}{last_mod} cmp $pipelines{$a}{last_mod} } keys %pipelines ) {
        push @pipelines_sorted, $pipelines{$pipeline};
        $pipelines_sorted[-1]{last_mod} = localtime( $pipelines_sorted[-1]{last_mod} );
    }
}

## populate the template
&print_template();

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

sub print_template {
    my $submenu_links = [
                            { label => 'new pipeline', is_last => 0, url => "./build_pipeline.cgi?repository_root=$repository_root" },
                        ];

    if ( $view eq 'component' ) {
        $tmpl->param( COMPONENTS        => \@components_sorted );
        push @$submenu_links, { label => 'view by pipeline', is_last => 0, url => "./pipeline_list.cgi?repository_root=$repository_root" };
        push @$submenu_links, { label => 'view by group', is_last => 1, url => "./pipeline_list.cgi?repository_root=$repository_root&view=group" };
    } elsif ( $view eq 'group' ) {
        $tmpl->param( PIPELINE_GROUPS   => \@pipeline_groups_sorted );
        push @$submenu_links, { label => 'view by pipeline', is_last => 0, url => "./pipeline_list.cgi?repository_root=$repository_root" };
        push @$submenu_links, { label => 'view by component', is_last => 0, url => "./pipeline_list.cgi?repository_root=$repository_root&view=component" };
        push @$submenu_links, { label => 'edit groups', is_last => 1, url => "javascript:showGroupModificationMenu()" };
    } else {
        $tmpl->param( PIPELINES        => \@pipelines_sorted );
        push @$submenu_links, { label => 'view by component', is_last => 0, url => "./pipeline_list.cgi?repository_root=$repository_root&view=component" };
        push @$submenu_links, { label => 'view by group', is_last => 1, url => "./pipeline_list.cgi?repository_root=$repository_root&view=group" };
    }

    ## populate the template with the values that will always be passed.
    $tmpl->param( PROJECT_CODE     => $project_code );
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
    $tmpl->param( PAGE_LINKS       => $page_links );
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

sub record_error {
    my $error_string = shift;
    
    $errors_found++;
    push @$error_msgs, { msg => $error_string };
}
