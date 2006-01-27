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

my $tmpl = HTML::Template->new( filename => 'templates/index.tmpl',
                                die_on_bad_params => 1,
                              );

## read the ergatis config file
my $ergatis_cfg = new Ergatis::ConfigFile( -file => "ergatis.ini" );
my $temp_space  = $ergatis_cfg->val( 'paths', 'temp_space' ) || die "temp_space not defined in ergatis.ini file";

my $cfg_md5 = $ergatis_cfg->project_list_md5();

## where should we find (or save) the running pipe list?
#my $running_pipeline_storable_path = "$temp_space/


## build the project list
my $registered_projects = [];
for my $label ( sort $ergatis_cfg->Parameters('projects') ) {
    push @$registered_projects, { 
                                    label => $label,
                                    repository_root => $ergatis_cfg->val('projects', $label),
                                };
}

my $running_pipelines = [];

## is there a stored version of the pipeline list?
my $running_dump_file = "$temp_space/$cfg_md5.ergatis.running.dump";

## if it exists and is less than 10 minutes old just display it.
if ( -e $running_dump_file && -M $running_dump_file < 0.0069444444) {
    $running_pipelines = retrieve $running_dump_file;

} else {
    &get_currently_running_pipelines( $running_pipelines );
    store $running_pipelines, $running_dump_file;
}


$tmpl->param( REGISTERED_PROJECTS => $registered_projects );
$tmpl->param( RUNNING_PIPELINES   => $running_pipelines );

print $tmpl->output;

sub get_currently_running_pipelines {
    my $list = shift;
    
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
            my $start_time = 'unknown';
            my $end_time = 'unknown';
            my $run_time = 'unknown';
            my $pipeline_user = 'unknown';
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

            next unless ( $state eq 'running' );

            ($start_time, $end_time, $run_time) = &time_info( $commandSet );
            my $filestat = stat($pipeline_file);
            $pipeline_user = getpwuid($filestat->uid);
            $last_mod = $filestat->mtime;

            ## this is done as a new twig parse since elements can be nested
            ## at any level.
            my %components = &component_count_hash( $pipeline_file );

            push @$list, {
                            label           => $label,
                            pipeline_id     => $pipeline_id,
                            state           => $state,
                            last_mod        => $last_mod,
                            run_time        => $run_time,
                            pipeline_user   => $pipeline_user,
    #                        components      => \@$component_aref,
    #                        component_count => $component_count,
    #                        component_label => $component_label,
    #                        view_link       => $view_link,
    #                        edit_link       => $edit_link,
                };
        }
    }    
}
