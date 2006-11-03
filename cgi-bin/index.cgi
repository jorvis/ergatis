#!/usr/local/bin/perl -w

use strict;
use CGI;
use CGI::Carp qw(fatalsToBrowser);
use Data::Dumper;
use Ergatis::Common;
use Ergatis::ConfigFile;
use HTML::Template;
use Monitor;
use Storable;
use XML::Twig;

my $q = new CGI;
print $q->header( -type => 'text/html' );

umask(0000);

## this toggle will force a rescan of the pipelines rather than pulling from 
##  the storable object.
my $update_cache = $q->param('update_cache') || 0;

my $tmpl = HTML::Template->new( filename => 'templates/index.tmpl',
                                die_on_bad_params => 0,
                              );

## read the ergatis config file
my $ergatis_cfg = new Ergatis::ConfigFile( -file => "ergatis.ini" );

my $temp_space = $ergatis_cfg->val( 'paths', 'temp_space' ) || die "temp_space not defined in ergatis.ini file";

## an MD5 is used on the project list so that many installations can
## share the same storable object.
my $cfg_md5 = $ergatis_cfg->project_list_md5();

my $table_cache_update_time = $ergatis_cfg->val( 'display_settings', 'pipeline_list_cache_time' ) || 10;  # in minutes
   $table_cache_update_time /= ( 60 * 24 );  ## convert to number of days (usually a decimal)

my $active_pipeline_age = $ergatis_cfg->val( 'display_settings', 'active_pipeline_age') || 24;

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
$tmpl->param( QUICK_LINKS         => &get_quick_links($ergatis_cfg) );
$tmpl->param( SUBMENU_LINKS       => [
                                        { label => 'update cache', is_last => 1, url => './index.cgi?update_cache=1' },
                                     ] );

print $tmpl->output;

sub get_pipeline_lists {
    my $running_list = [];
    my $active_list  = [];
    
    for ( @$registered_projects ) {
        my $label           = $_->{label};
        my $repository_root = $_->{repository_root};

        ## CATCH WARNING HERE LATER
        next unless ( -d "$repository_root/workflow/runtime/pipeline" );

        ## open the pipeline dir
        opendir ( my $idh, "$repository_root/workflow/runtime/pipeline" ) || die "can't read directory $repository_root/workflow/runtime/pipeline: $!";

        for my $pipeline_id ( readdir $idh ) {
            my $state = '';
            my $last_mod = '';
            my $pipeline_file = "$repository_root/workflow/runtime/pipeline/$pipeline_id/pipeline.xml.instance";
            
            ## it may have been compressed
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
            
            my $pipeline_file_fh;
            if ($pipeline_file =~ /\.gz/) {
                open($pipeline_file_fh, "<:gzip", "$pipeline_file") || die "can't read $pipeline_file: $!"; 
            } else {
                open($pipeline_file_fh, "<$pipeline_file") || die "can't read $pipeline_file: $!";       
            }

            my $twig = new XML::Twig;
            
            ## this would otherwise die if the pipeline file is 0 bytes
            if (! eval { $twig->parse($pipeline_file_fh) } ) {
                ## note, this would be better handled with some message showing that the file empty
                next;
            }

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

            ## depending on the state, grab the top-level error message if there is one
            ##  there's no use doing this for running pipelines, since workflow won't
            ##  propagate the error until it's finished.
            my $error_message = 0;
            if ( $state eq 'error' || $state eq 'failed' ) {
                if ( $commandSet->has_child('status') && 
                     $commandSet->first_child('status')->has_child('message') ) {

                    $error_message = $commandSet->first_child('status')->first_child('message')->text();
                    
                    ## handle illegal characters
                    $error_message = CGI::escapeHTML($error_message);
                }
            }
            
            my %components = &component_count_hash( $pipeline_file );
            my $component_count = 0;
            my $component_aref;
            foreach my $component (sort keys %components) {
                $component_count += $components{$component}{count};
                push @$component_aref, { name => $component, 
                                         count => $components{$component}{count},
                                         error_count => $components{$component}{error_count},
                                       };
            }
            
            ## reformat component_count to include a label
            my $component_label = ' component';
            if ($component_count != 1) {
                $component_label = ' components';
            }
            
            my $last_mod_secs = $last_mod;
            $last_mod = localtime($last_mod);
            
            my $pipeline_info = {
                            label           => $label,
                            project_url     => "./pipeline_list.cgi?repository_root=$repository_root",
                            pipeline_id     => $pipeline_id,
                            pipeline_url    => "./view_pipeline.cgi?instance=$pipeline_file",
                            state           => $state,
                            last_mod        => $last_mod,
                            last_mod_secs   => $last_mod_secs,
                            run_time        => $run_time,
                            pipeline_user   => $pipeline_user,
                            components      => \@$component_aref,
                            component_count => $component_count,
                            component_label => $component_label,
                            error_message   => $error_message,
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
    my @sortedrunning_list = sort {$b->{'last_mod_secs'} cmp $a->{'last_mod_secs'}} @$running_list;    
    my @sortedactive_list = sort {$b->{'last_mod_secs'} cmp $a->{'last_mod_secs'}} @$active_list;
    return ( \@sortedrunning_list, \@sortedactive_list );
}
