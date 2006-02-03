#!/usr/local/bin/perl -w

use strict;
use CGI;
use CGI::Carp qw(fatalsToBrowser);
use Data::Dumper;
use Ergatis::ConfigFile;
use HTML::Template;
use Monitor;
use Storable;
use XML::Twig;

my $q = new CGI;
print $q->header( -type => 'text/html' );

## this toggle will force a rescan of the pipelines rather than pulling from 
##  the storable object.
my $update_cache = $q->param('update_cache') || 0;

my $tmpl = HTML::Template->new( filename => 'templates/index.tmpl',
                                die_on_bad_params => 1,
                              );

## read the ergatis config file
my $ergatis_cfg = new Ergatis::ConfigFile( -file => "ergatis.ini" );

my $temp_space          = $ergatis_cfg->val( 'paths', 'temp_space' ) || die "temp_space not defined in ergatis.ini file";

## an MD5 is used on the project list so that many installations can
## share the same storable object.
my $cfg_md5 = $ergatis_cfg->project_list_md5();

my $table_cache_update_time = $ergatis_cfg->val( 'settings', 'pipeline_list_cache_time' ) || 10;  # in minutes
   $table_cache_update_time /= ( 60 * 24 );  ## convert to number of days (usually a decimal)

my $active_pipeline_age = $ergatis_cfg->val( 'settings', 'active_pipeline_age') || 24;

## build the project list
my $registered_projects = [];
for my $label ( sort $ergatis_cfg->Parameters('projects') ) {
    push @$registered_projects, { 
                                    label => $label,
                                    repository_root => $ergatis_cfg->val('projects', $label),
                                };
}

my $running_pipelines = [];
my $active_pipelines = [];

## is there a stored version of the running pipeline list?
my $running_dump_file = "$temp_space/$cfg_md5.ergatis.running.dump";
my $active_dump_file  = "$temp_space/$cfg_md5.ergatis.active.dump";

## if it exists and is less than the cache update time just display it (unless we're forcing an update).
my $cache_file_age = -M $running_dump_file;

if ( -e $running_dump_file && 
     -e $active_dump_file && 
     $cache_file_age < $table_cache_update_time && 
     ! $update_cache ) {
     
    $running_pipelines = retrieve $running_dump_file;
    $active_pipelines  = retrieve $active_dump_file;

} else {
    ( $running_pipelines, $active_pipelines ) = get_pipeline_lists();
    store $running_pipelines, $running_dump_file;
    store $active_pipelines, $active_dump_file;
    $cache_file_age = -M $running_dump_file;
}


$tmpl->param( REGISTERED_PROJECTS => $registered_projects );
$tmpl->param( RUNNING_PIPELINES   => $running_pipelines );
$tmpl->param( ACTIVE_PIPELINES    => $active_pipelines );
$tmpl->param( ACTIVE_PIPELINE_AGE => $active_pipeline_age );
$tmpl->param( CACHE_FILE_AGE      => int($cache_file_age * 1440) );

print $tmpl->output;

sub get_pipeline_lists {
    my $running_list = [];
    my $active_list  = [];
    
    for ( @$registered_projects ) {
        my $label           = $_->{label};
        my $repository_root = $_->{repository_root};

        ## CATCH WARNING HERE LATER
        next unless ( -d "$repository_root/Workflow/pipeline" );

        ## open the pipeline dir
        opendir ( my $idh, "$repository_root/Workflow/pipeline" ) || die "can't read directory $repository_root/Workflow/pipeline: $!";

        for my $pipeline_id ( readdir $idh ) {
            my $state = '';
            my $last_mod = '';
            my $pipeline_file = "$repository_root/Workflow/pipeline/$pipeline_id/pipeline.xml.instance";

            next unless ( -e $pipeline_file );

            my $twig = new XML::Twig;
            $twig->parsefile($pipeline_file);

            my $commandSetRoot = $twig->root;
            my $commandSet = $commandSetRoot->first_child('commandSet');

            next if (! $commandSet );

            if ( $commandSet->first_child('state') ) {
                $state  = $commandSet->first_child('state')->text();
            }

            my $filestat = stat($pipeline_file);

            $last_mod = $filestat->mtime;
            
            ## check the time here.  we'll skip this one unless
            ## it is either running or less than active_pipeline_age time
            next unless ( $state eq 'running' || (time - $last_mod) < ($active_pipeline_age * 3600) );
            
            
            my $pipeline_user = getpwuid($filestat->uid);
            my ($start_time, $end_time, $run_time) = &time_info( $commandSet );
            
            my %components = &component_count_hash( $pipeline_file );
            my $component_count = 0;
            my $component_aref;
            foreach my $component (sort keys %components) {
                $component_count += $components{$component};
                push @$component_aref, { name => $component, count => $components{$component} };
            }
            
            ## reformat component_count to include a label
            my $component_label = ' component';
            if ($component_count != 1) {
                $component_label = ' components';
            }
            
            $last_mod = localtime($last_mod);
            
            my $pipeline_info = {
                            label           => $label,
                            project_url     => "./pipeline_list.cgi?repository_root=$repository_root",
                            pipeline_id     => $pipeline_id,
                            pipeline_url    => "./view_workflow_pipeline.cgi?instance=$pipeline_file",
                            state           => $state,
                            last_mod        => $last_mod,
                            run_time        => $run_time,
                            pipeline_user   => $pipeline_user,
                            components      => \@$component_aref,
                            component_count => $component_count,
                            component_label => $component_label,
                };
                
            ## if the pipeline is still running, add it to that list.  else, add it to the 'active'
            ##   list (since we filtered previously)
            if ( $state eq 'running' ) {
                push @$running_list, $pipeline_info;
            } else {
                push @$active_list, $pipeline_info;
            }
        }
    }    
    
    return ( $running_list, $active_list );
}
