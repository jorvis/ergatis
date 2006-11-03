#!/usr/local/bin/perl -w

use strict;
use CGI;
use CGI::Carp qw(fatalsToBrowser);
use Ergatis::Common;
use Ergatis::ConfigFile;
use HTML::Template;
use Monitor;
use XML::Twig;

my $q = new CGI;

print $q->header( -type => 'text/html' );
my $tmpl = HTML::Template->new( filename => 'templates/pipeline_list.tmpl',
                                die_on_bad_params => 1,
                              );

my $repository_root = $q->param("repository_root") || die "pass a repository root";

## quota information (only works if in /usr/local/annotation/SOMETHING)
my $quotastring = &quota_string($repository_root);

my $ergatis_cfg = new Ergatis::ConfigFile( -file => "ergatis.ini" );
my $display_codebase = $ergatis_cfg->val( 'display_settings', 'display_codebase') || 0;

my $pipeline_root = "$repository_root/workflow/runtime/pipeline";
my $shared_conf_path = "$repository_root/workflow/project.config";
my $errors_found  = 0;
my $error_msgs = [];
my %pipelines;
my $pipeline_count = 0;

## make sure the directory exists
if (! -e $pipeline_root) {
    &record_error("$pipeline_root not found.  enter a repository root that contains a workflow directory.");
    &print_template();
}

my $rdh;
if (! opendir $rdh, "$pipeline_root" ) {
    &record_error("can't open pipeline root directory: $!");
    &print_template();
}

## make sure it has a 'workflow' directory
if (! -d "$repository_root/workflow" ) {
    &record_error("directory $repository_root/workflow not found");
}

## make sure it has a 'workflow/lock_files' directory
if (! -d "$repository_root/workflow/lock_files" ) {
    &record_error("directory $repository_root/workflow/lock_files not found");

} else {
    ## make sure it is writeable
    if (! -w "$repository_root/workflow/lock_files") {
        &record_error("$repository_root/workflow/lock_files not writable");
    }
}

## make sure it has a 'workflow/project_id_repository' directory
if (! -d "$repository_root/workflow/project_id_repository" ) {
    &record_error("directory $repository_root/workflow/project_id_repository not found");
} else {
    ## make sure it is writeable
    if (! -w "$repository_root/workflow/project_id_repository") {
        &record_error("$repository_root/workflow/project_id_repository not writable");
    }
}

## make sure it has a shared conf file
if (! -e $shared_conf_path ) {
    &record_error("$shared_conf_path not found");
}

## quit and throw the template if there were errors.
if ( scalar @{$error_msgs} ) {
    &print_template();
}

## pull the ergatis dir from the shared conf file
my $shared_conf = new Ergatis::ConfigFile( -file => $shared_conf_path );
my $ergatis_dir = $shared_conf->val('project', '$;ERGATIS_DIR$;') || 'unknown';

foreach my $pipeline_id ( readdir $rdh ) {
    next unless ( $pipeline_id =~ /^\d+$/ );
    
    print STDERR "parsing pipeline id $pipeline_id\n";
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
    my $pipeline_file = "$repository_root/workflow/runtime/pipeline/$pipeline_id/pipeline.xml";  ## may be modified below
    my $archive_link = "./archive_pipeline_form.cgi?repository_root=$repository_root&amp;pipeline_id=$pipeline_id";
    my $links_enabled = 1;
    my $error_message = 0;
    
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

        ## if only the pipeline.xml exists, we can do less
        if (! -e $pipeline_file ) {

            if ( -e "$pipeline_file.gz" ) {
                $pipeline_file .= '.gz';
            } else {
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
        my %components = &component_count_hash( $pipeline_file );

        foreach my $component (sort keys %components) {
            $component_count += $components{$component}{count};
            push @$component_aref, { name => $component, 
                                     count => $components{$component}{count},
                                     error_count => $components{$component}{error_count},
                                   };
        }

        ## reformat component_count to include a label
        if ($component_count == 1) {
            $component_label = ' component';
        }
    }
    
    my $view_link = "./view_pipeline.cgi?instance=$pipeline_file";
    my $edit_link = "./show_pipeline.cgi?xmltemplate=$pipeline_file&amp;edit=1";
    
    $pipelines{$pipeline_id} = { 
                        pipeline_id     => $pipeline_id,
                        state           => $state,
                        last_mod        => $last_mod,
                        run_time        => $run_time,
                        pipeline_user   => $pipeline_user,
                        components      => \@$component_aref,
                        component_count => $component_count,
                        component_label => $component_label,
                        view_link       => $view_link,
                        edit_link       => $edit_link,
                        archive_link    => $archive_link,
                        links_enabled   => $links_enabled,
                        error_message   => $error_message,
                      };
}

## sort the pipelines
my @pipelines_sorted;

for my $pipeline ( sort { $pipelines{$b}{last_mod} cmp $pipelines{$a}{last_mod} } keys %pipelines ) {
    $pipelines{$pipeline}{last_mod} = localtime( $pipelines{$pipeline}{last_mod} );
    push @pipelines_sorted, $pipelines{$pipeline};
}

## populate the template
&print_template();

sub print_template {

    ## populate the template with the values that will always be passed.
    $tmpl->param( ERRORS_FOUND     => $errors_found );
    $tmpl->param( REPOSITORY_ROOT  => $repository_root );
    $tmpl->param( ERROR_MSGS       => $error_msgs );
    $tmpl->param( PIPELINES        => \@pipelines_sorted );
    $tmpl->param( PIPELINE_COUNT   => $pipeline_count );
    $tmpl->param( QUOTA_STRING     => $quotastring );
    $tmpl->param( ERGATIS_DIR      => $ergatis_dir );
    $tmpl->param( DISPLAY_CODEBASE => $display_codebase );
    $tmpl->param( QUICK_LINKS         => &get_quick_links($ergatis_cfg) );
    $tmpl->param( SUBMENU_LINKS       => [
                                            { label => 'new pipeline', is_last => 1, url => "./build_pipeline.cgi?repository_root=$repository_root" },
                                         ] );

    ## print the template
    print $tmpl->output();
    exit;
}

sub record_error {
    my $error_string = shift;
    
    $errors_found++;
    push @$error_msgs, { msg => $error_string };
}
