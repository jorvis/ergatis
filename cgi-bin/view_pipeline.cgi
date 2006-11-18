#!/usr/local/bin/perl -w

use strict;
use CGI;
use CGI::Carp qw(fatalsToBrowser);
use Date::Manip;
use Ergatis::Common;
use Ergatis::ConfigFile;
use Ergatis::Monitor;
use HTML::Template;
use XML::Twig;


my $q = new CGI;
print $q->header( -type => 'text/html' );

my $tmpl = HTML::Template->new( filename => 'templates/view_pipeline.tmpl',
                                die_on_bad_params => 1,
                              );

## read the ergatis config file
my $ergatis_cfg = new Ergatis::ConfigFile( -file => "ergatis.ini" );

my $repository_root;
my $xml_input = $q->param("instance") || die "pass instance";

my $project = '?';
my $pipelineid;

## make sure the file exists, check for .gz version
if (! -f $xml_input ) {
    if (-f "$xml_input.gz") {
        $xml_input .= '.gz';
    } else {
        die "$xml_input doesn't exist";
    }
}

## extract the project
if ( $xml_input =~ m|(.+/(.+?))/workflow/runtime/pipeline/(\d+)/| ) {
    $repository_root = $1;
    $project = $2;
    $pipelineid = $3;
} else {
    die "failed to extract a repository_root from $xml_input.  expected a workflow/runtime subdirectory somewhere."
}

my $file = $xml_input;

my $ifh;
if ($file =~ /\.gz/) {
    open($ifh, "<:gzip", "$file") || die "can't read $file: $!"; 
} else {
    open($ifh, "<$file") || die "can't read $file: $!";       
}

my $twig = new XML::Twig;

$twig->parse($ifh);

my $commandSetRoot = $twig->root;
my $commandSet = $commandSetRoot->first_child('commandSet');

## pull desired info out of the root commmandSet
## pull: project space usage
my ($starttime, $endtime, $lastmodtime, $state, $runtime) = ('n/a', 'n/a', '', 'unknown', 'n/a');

if ( $commandSet->first_child('state') ) {
    $state  = $commandSet->first_child('state')->text();
}

($starttime, $endtime, $runtime) = &time_info($commandSet);

my $filestat = stat($file);
my $user = getpwuid($filestat->uid);
$lastmodtime = time - $filestat->mtime;
$lastmodtime = strftime( "%H hr %M min %S sec", reverse split(/:/, DateCalc("today", ParseDate(DateCalc("now", "- ${lastmodtime} seconds")) ) ));

## quota information (only works if in /usr/local/annotation/SOMETHING)
my $quotastring = &quota_string($repository_root);

## if state is not complete, start an update timer
if ($state ne 'complete') {
   #print "<script>pipelineCountdown();</script>\n";
}

## i'm just building a string here for the component html string because
##  i don't see how to do it with HTML::Template since it doesn't support
##  any kind of recursion.  (since components can be nested an any 
##  arbitrary level within command sets)
my $component_html = '';

## parse the root commandSet to find components and set definitions
parseCommandSet( $commandSet );


$tmpl->param( PIPELINE_FILE       => $file );
$tmpl->param( START_TIME          => $starttime );
$tmpl->param( END_TIME            => $endtime );
$tmpl->param( LAST_MOD_TIME       => $lastmodtime );
$tmpl->param( PIPELINE_STATE      => $state );
$tmpl->param( USER                => $user );
$tmpl->param( RUNTIME             => $runtime );
$tmpl->param( PROJECT             => $project );
$tmpl->param( QUOTA_STRING        => $quotastring );
$tmpl->param( PIPELINE_ID         => $pipelineid );
$tmpl->param( COMPONENT_HTML      => $component_html );
$tmpl->param( PAGE_TITLE          => "$project|$state|pipeline" );
$tmpl->param( QUICK_LINKS         => &get_quick_links($ergatis_cfg) );
$tmpl->param( SUBMENU_LINKS       => [
                                        { label => 'pipeline list', is_last => 0, url => "./pipeline_list.cgi?repository_root=$repository_root" },
                                        { label => 'new pipeline', is_last => 0, url => "./build_pipeline.cgi?repository_root=$repository_root" },
                                        { label => 'rerun', is_last => 0, url => "./run_pipeline.cgi?pipeline_xml=$file&pipeline_id=$pipelineid&rerun=1" },
                                        { label => 'kill', is_last => 0, url => "./kill_wf.cgi?instancexml=$file" },
                                        { label => 'view xml', is_last => 1, url => "./view_formatted_xml_source.cgi?file=$file" },
                                     ] );

print $tmpl->output;

exit(0);



sub parseCommandSet {
    my ($commandSet) = @_;
   
    my $commandSet_name = $commandSet->first_child('name')->text();
    my $type = $commandSet->{att}->{type};

    ## there are three cases here.  it's either the start, a set, or a component
    ## check for the pipeline start
    if ( $commandSet_name =~ /^start pipeline:(.*)/ ) {
        parseCommandSetChildren( $commandSet );
    
    ## serial or parallel sets will be named 'serial' or 'parallel'
    } elsif ( $commandSet_name eq 'serial' || $commandSet_name eq 'parallel' ) {

        $component_html .= "<ul class='$type'>\n";

        ## some special elements are needed for parallel sets
        if ($type eq 'parallel') {
            $component_html .= "    <h1 class='commandset_type'>parallel group</h1>\n";
        } else {
            $component_html .= "    <h1 class='commandset_type'>serial group</h1>\n";
        }

        ## look at each child.
        parseCommandSetChildren( $commandSet );
        $component_html .= "</ul>\n";

    } elsif ( $commandSet_name =~ /^((.+?)\.(.+))/) {
        my $filebased_subflow = '';
        my $name_token = $1;

        ## this should be changed to 0 if the interface shouldn't auto-update the 
        #   box for this component.
        my $do_auto_update = 1;
        my $user_msg = 'state: wait for update';
        
        ## this is a component, so the XML here will consist of the following:
        #   command: replace_config_keys
        #   command: replace_template_keys
        #   commandSet: has a file-based subflow that is the component.xml

        ## TODO: check the replace_config_keys step
        ## TODO: check the replace_template_keys step
        
        my $subcommandSet = $commandSet->first_child('commandSet') || 0;
        if ( $subcommandSet ) {
            my $fileName = $subcommandSet->first_child('fileName') || 0;
            if ($fileName) {
                $filebased_subflow = $fileName->text;
            } else {
                ## we expected a fileName here
                ##  TODO: handle error?
            }
        } else {
            ## we expected a commandSet here.  
            ##  TODO: handle error?
        }
        
        $component_html .= <<ComponeNTBlock;
<ul class='component' id='$name_token'>
    <li>
        <div class="component_label"><b>component:</b> $name_token</div>
    </li>
    <li>state: $user_msg</li>
</ul>
ComponeNTBlock

        ## should we spawn an autoupdate?
        if ( $do_auto_update ) {
            $component_html .= "<script>sendComponentUpdateRequest('./component_summary.cgi?pipeline=$filebased_subflow&ul_id=$name_token&parent_pipeline=$xml_input', updateComponent, '$name_token', '$filebased_subflow', '$xml_input');</script>\n";
        }
    }
}


sub parseCommandSetChildren {
    my $commandSet = shift;
    
    ## look at each child.
    for my $child ( $commandSet->children() ) {
        ## if it is a command set, parse through it. 
        if ($child->gi eq 'commandSet') {
            parseCommandSet($child);
        }
    }
}


exit;
