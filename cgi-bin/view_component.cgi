#!/usr/local/bin/perl -w

use strict;
use CGI;
use CGI::Carp qw(fatalsToBrowser);
use Date::Manip;
use Ergatis::Common;
use Ergatis::ConfigFile;
use File::stat;
use HTML::Template;
use Monitor;
use POSIX;
use XML::Twig;


my $q = new CGI;

print $q->header( -type => 'text/html' );

my $tmpl = HTML::Template->new( filename => 'templates/view_component.tmpl',
                                die_on_bad_params => 1,
                              );

## read the ergatis config file
my $ergatis_cfg = new Ergatis::ConfigFile( -file => "ergatis.ini" );


my $pipeline_xml = $q->param("pipeline_xml") || die "pass pipeline_xml";

## it may have been compressed
if (! -e $pipeline_xml && -e "$pipeline_xml.gz") {
    $pipeline_xml .= '.gz';
}

my $pipeline_xml_fh;
if ($pipeline_xml =~ /\.gz/) {
    open($pipeline_xml_fh, "<:gzip", "$pipeline_xml") || die "can't read $pipeline_xml: $!"; 
} else {
    open($pipeline_xml_fh, "<$pipeline_xml") || die "can't read $pipeline_xml: $!";       
}

my $twig = XML::Twig->new( );
$twig->parse($pipeline_xml_fh);

my $parent_commandset = $twig->root->first_child('commandSet');
my $component_state = $parent_commandset->first_child('state')->text || 'unknown';

my $component_name = 'component';
my $project = '?';
my $repository_root = '';
my $pipeline_id = '';
my $output_token = '';

if ( $pipeline_xml =~ m|(.+/(.+?))/workflow/runtime/(.+?)/(\d+)_(.+?)/| ) {
    $repository_root = $1;
    $project = $2;
    $component_name = $3;
    $pipeline_id = $4;
    $output_token = $5;
}

my $parent_pipeline = '';
if ( $parent_commandset->first_child('parentFileName') ) {
    $parent_pipeline = $parent_commandset->first_child('parentFileName')->text;
}

my %lockfile_info = &parse_pipeline_run_lock_file( "$repository_root/workflow/lock_files/pid.$pipeline_id" );

## quota information (only works if in /usr/local/annotation/SOMETHING)
my $quotastring = &quota_string($repository_root);

my ($starttime, $endtime, $lastmodtime, $state, $runtime) = ('n/a', 'n/a', '', 'unknown', 'n/a');

($starttime, $endtime, $runtime) = &time_info($parent_commandset);

my $filestat = stat($pipeline_xml);
my $user = getpwuid($filestat->uid);
$lastmodtime = time - $filestat->mtime;
$lastmodtime = strftime( "%H hr %M min %S sec", reverse split(/:/, DateCalc("today", ParseDate(DateCalc("now", "- ${lastmodtime} seconds")) ) ));

my $component_html = "<div id='workflowcontainer'>\n";

my $elements = [];

## look at each of the children of the root
foreach my $child ( $parent_commandset->children() ) {
    
    ## if this is a command, process its details
    if ($child->gi eq 'command') {
        
        my %parts = &process_command($twig, $child);
        push @{$elements}, \%parts;
    
    ## if it is a commandSet, it should be a file-based subflow
    } elsif ($child->gi eq 'commandSet') {
    
        ## make sure it has a fileName element
        if ( $child->has_child('fileName') ) {
        
            push @$elements, { is_command => 0 };
        
            ## set a label for this
            $$elements[-1]->{subflow_label} = $child->first_child('name')->text();
            
            &parse_iterator_xml( $child->first_child('fileName')->text() );
        }
    }
}

$tmpl->param( PIPELINE_FILE       => $pipeline_xml );
$tmpl->param( START_TIME          => $starttime );
$tmpl->param( END_TIME            => $endtime );
$tmpl->param( LAST_MOD_TIME       => $lastmodtime );
$tmpl->param( PIPELINE_STATE      => $state );
$tmpl->param( USER                => $user );
$tmpl->param( RUNTIME             => $runtime );
$tmpl->param( RETRIES             => $lockfile_info{retries} );
$tmpl->param( EXEC_USER           => $lockfile_info{execuser} );
$tmpl->param( HOST_NAME           => $lockfile_info{hostname} );
$tmpl->param( PID                 => $lockfile_info{pid} );
$tmpl->param( PROJECT             => $project );
$tmpl->param( QUOTA_STRING        => $quotastring );
$tmpl->param( PIPELINE_ID         => $pipeline_id );
#$tmpl->param( COMPONENT_HTML      => $component_html );

$tmpl->param( ELEMENTS            => $elements );

$tmpl->param( PAGE_TITLE          => "$project|$component_name|$state" );
$tmpl->param( QUICK_LINKS         => &get_quick_links($ergatis_cfg) );
$tmpl->param( SUBMENU_LINKS       => [
                                        { label => 'pipeline view', is_last => 0, url => "./view_pipeline.cgi?instance=$parent_pipeline" },
                                        { label => 'view configuration', is_last => 0, url => 
                                        "./view_formatted_ini_source.cgi?file=$repository_root/workflow/runtime/$component_name/${pipeline_id}_$output_token/$component_name.$output_token.final.config" },
                                        { label => 'view xml', is_last => 1, url => "./view_formatted_xml_source.cgi?file=$pipeline_xml" },
                                     ] );

print $tmpl->output;

exit(0);



sub parse_iterator_xml {
    my $filename = shift;
    
    ## it may have been compressed
    if (! -e $filename) {
        if (-e "$filename.gz") {
            $filename .= '.gz';
        } else {
            $component_html .= "            <div>not yet created</div>\n";    
        }
    }
    
    ## make sure this is a iN.xml file
    if ( $filename !~ /i\d+.xml/ ) {
        $component_html .= "            <div>unable to handle $filename</div>\n";
        return;
    }
    
    my ($component_start_time, $component_end_time, $component_state);
    
    my $filename_fh;
    if ($filename =~ /\.gz/) {
        open($filename_fh, "<:gzip", "$filename") || die "can't read $filename: $!"; 
    } else {
        open($filename_fh, "<$filename") || die "can't read $filename: $!";       
    }
    
    ## create the twig
    my $twig = XML::Twig->new( twig_roots => {
                                    'commandSet' => \&process_subflowgroup,
                               }
                             );
    $twig->parse($filename_fh);
    
}

sub process_subflowgroup {
    my ($twig, $commandSet) = @_;
   
    ## within the iN.xml file, each commandSet with a fileName defines a gN.xml file.
    
    ## do nothing unless this commandSet has a file-based subflow
    return unless $commandSet->has_child('fileName');
    
    my %sg_props = (
                      execution_host => '',
                      file => '',
                      grid_id => '?',
                      group_num => '?',
                      message => '',
                      name => '',
                      ret_value => 'unknown',
                      state => 'unknown',
                      workflow_id => '?',
                   );
    
    ## get the pointer to the gN file
    $sg_props{file} = $commandSet->first_child('fileName')->text();
    
    ## pull the group number out of the file name:
    if ( $sg_props{file} =~ m|i(\d+)/g\d+/g(\d+).xml| ) {
        $sg_props{name} = "$1$2";
        $sg_props{group_num} = $2;
    }
    
    ## get the state, if it has one
    if ( $commandSet->first_child('state') ) {
        $sg_props{state} = $commandSet->first_child('state')->text();
    }
    
    ## grab data from the dceSpec if it has one
    if ( $commandSet->first_child('dceSpec') ) {
        if ( $commandSet->first_child('dceSpec')->first_child('executionHost') ) {
            $sg_props{execution_host} = $commandSet->first_child('dceSpec')->first_child('executionHost')->text();
        }
        
        if ( $commandSet->first_child('dceSpec')->first_child('jobID') ) {
            $sg_props{grid_id} = $commandSet->first_child('dceSpec')->first_child('jobID')->text();
        }
    }
    
    if ( $commandSet->first_child('id') ) {
        $sg_props{workflow_id} = $commandSet->first_child('id')->text();
    }
    
    ( $sg_props{start_time}, $sg_props{end_time}, $sg_props{run_time} ) = &time_info($commandSet);
    
    ## if there is a status and a message, grab it
    if ( $commandSet->first_child('status') ) {
        if ( $commandSet->first_child('status')->first_child('retValue') ) {
            $sg_props{ret_value} = $commandSet->first_child('status')->first_child('retValue')->text;
        }
        
        if ( $commandSet->first_child('status')->first_child('message') ) {
            $sg_props{message} = $commandSet->first_child('status')->first_child('message')->text;
        }
        
        ## don't include 'command set ... finished' messages
        $sg_props{message} =~ s/command set .* finished//i;
    }

    $sg_props{hostsrvstr} = join(',',split(/\./,$sg_props{execution_host}));

    push @{$$elements[-1]->{sg_props}}, \%sg_props;
}

