#!/usr/bin/perl -w

use strict;
use CGI;
use CGI::Carp qw(fatalsToBrowser);
use Ergatis::Common;
use Ergatis::ConfigFile;
use Ergatis::Monitor;
use HTML::Template;
use XML::Twig;

my $q = new CGI;

print $q->header( -type => 'text/html' );
my $tmpl = HTML::Template->new( filename => 'templates/pipeline_list.tmpl',
                                die_on_bad_params => 0,
                                global_vars => 1
                              );

my $repository_root = $q->param("repository_root") || die "pass a repository root";

## quota information (only works if in /usr/local/annotation/SOMETHING)
my $quotastring = &quota_string($repository_root);

my $ergatis_cfg = new Ergatis::ConfigFile( -file => "ergatis.ini" );
my $display_codebase = $ergatis_cfg->val( 'display_settings', 'display_codebase') || 0;

my $pipeline_root = "$repository_root/workflow/runtime/pipeline";
my $project_conf_path = "$repository_root/workflow/project.config";
my $errors_found  = 0;
my $error_msgs = [];
my %pipelines;
my $pipeline_count = 0;
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

foreach my $pipeline_id ( readdir $rdh ) {
    next unless ( $pipeline_id =~ /^\d+$/ );
    
    $pipeline_count++;

    my $state = 'unknown';
    my $last_mod = 'unknown';
    my $start_time = 'unknown';
    my $end_time = 'unknown';
    my $run_time = 'unknown';
    my $pipeline_user = 'unknown';
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

        my $filestat = stat($pipeline_file);
        $pipeline_user = getpwuid($filestat->uid);
        $last_mod = $filestat->mtime;
        
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
                        pipeline_id     => $pipeline_id,
                        state           => $state,
                        is_running      => $state eq 'running' ? 1 : 0,
                        last_mod        => $last_mod,
                        run_time        => $run_time,
                        pipeline_user   => $pipeline_user,
                        components      => $component_aref,
                        component_count => $component_count,
                        component_label => $component_label,
                        view_link       => $view_link,
                        archive_link    => $archive_link,
                        links_enabled   => $links_enabled,
                        error_message   => $error_message,
                        has_comment     => length($pipeline_comment),
                        pipeline_comment => $pipeline_comment,
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
        }
    }    

    
    for my $group ( keys %pipeline_groups ) {
        push @pipeline_groups_sorted, { name => $group, pipelines => \@{$pipeline_groups{$group}} };
    }

} else {

    for my $pipeline ( sort { $pipelines{$b}{last_mod} cmp $pipelines{$a}{last_mod} } keys %pipelines ) {
        $pipelines{$pipeline}{last_mod} = localtime( $pipelines{$pipeline}{last_mod} );
        push @pipelines_sorted, $pipelines{$pipeline};
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
    $tmpl->param( PIPELINE_COUNT   => $pipeline_count );
    $tmpl->param( QUOTA_STRING     => $quotastring );
    $tmpl->param( ERGATIS_DIR      => $ergatis_dir );
    $tmpl->param( DISPLAY_CODEBASE => $display_codebase );
    $tmpl->param( COMPONENT_VIEW   => $view eq 'component' ? 1 : 0 );
    $tmpl->param( GROUP_VIEW       => $view eq 'group' ? 1 : 0 );
    $tmpl->param( QUICK_LINKS      => &get_quick_links($ergatis_cfg) );
    $tmpl->param( SUBMENU_LINKS    => $submenu_links );

    ## print the template
    print $tmpl->output();
    exit;
}

sub record_error {
    my $error_string = shift;
    
    $errors_found++;
    push @$error_msgs, { msg => $error_string };
}
