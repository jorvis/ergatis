#!/usr/bin/perl -w

use strict;
use CGI;
use CGI::Carp qw(fatalsToBrowser);
use Ergatis::ConfigFile;
use Ergatis::Monitor;
use HTML::Template;
use XML::Twig;
use File::Basename;

my $q = new CGI;

print $q->header( -type => 'text/html' );

my $tmpl = HTML::Template->new( filename => 'templates/component_summary.tmpl',
                                die_on_bad_params => 1,
                              );

## setting while testing.  if true, the parser will go down the the most finite command level
#   this takes too long on larger command sets, so i'm using it for testing now.
my $max_command_resolution = 0;

## will be like:
## /usr/local/annotation/TGA1/workflow/runtime/split_fasta/29134_test2/component.xml
my $pipeline = $q->param("pipeline") || die "pass pipeline";
my $p_pipeline_state = $q->param("parent_pipeline_state") || die "pass parent pipeline state";

## will be like:
## split_fasta.test2
my $ul_id = $q->param("ul_id") || die "pass ul_id";

## will be like:
## /usr/local/annotation/AA1/workflow/runtime/pipeline/29671/pipeline.xml
my $parent_pipeline = $q->param("parent_pipeline") || '';
my $pipeline_exists = 0;

my $progress_image_width = 500;
my $component_state = 'unknown';
my $command_count = 0;

my %states;
my $has_multiple_states = 0;
my @state_elements;

my @messages;
my %message_counts;
my $messages_line = '';
#time variables
my ($start_time, $end_time, $lastmodtime, $state, $runtime) = ('n/a', 'n/a', '', 'unknown', 0);
my $got_time_info = 0;
my $current_step = '';
my $update_interval = 61;

## we can parse some information out of the standardized instance file path
##  like: /usr/local/scratch/annotation/EHA1/workflow/runtime/iprscan/30835_default/component.xml
my ($root_dir, $project, $component, $token, $pipelineid, $component_conf_varreplaced, $component_conf_nonvarreplaced);
if ($pipeline =~ m|(.+/(.+)/workflow/runtime/(.+?)/([A-Z0-9]+)_(.+?))/component.xml|) {
    ($project, $component, $pipelineid, $token) = ($2, $3, $4, $5);
    $component_conf_nonvarreplaced = "$1/$component.$token.user.config";
    $component_conf_varreplaced = "$1/$component.$token.final.config";
} else {
    print "invalid instance file path format";
    exit(1);
}

if (-s $pipeline || -s "$pipeline.gz") {
    $pipeline_exists = 1;
    
    my $pipeline_fh;
    if ($pipeline =~ /\.gz/) {
        open($pipeline_fh, "<:gzip", "$pipeline") || die "can't read $pipeline: $!"; 
    } elsif ( -e "$pipeline.gz" ) {
        open($pipeline_fh, "<:gzip", "$pipeline.gz") || die "can't read $pipeline: $!"; 
    } else {
        open($pipeline_fh, "<$pipeline") || die "can't read $pipeline: $!";       
    }

    my $twig = new XML::Twig;
       $twig->parse($pipeline_fh);
    my $commandSetRoot = $twig->root;
    parseCommandSet( $commandSetRoot->first_child('commandSet'), $pipeline );
    
    ## build the line that lists each status and its count
    ## at the same time we can calculate the width of each state for the progress bar
    my $state_count = scalar keys %states;
    my $width_used = 0;  
    my $states_handled = 0;
   
    for my $status (sort keys %states) {
        ## each status gives a percentage of the total command_count
        my $width = int( ($states{$status} / $command_count) * $progress_image_width);
        $width_used += $width;
        $states_handled++;
        
        ## if this is the last state and there are unused pixels in the bar, just tag
        ## the unused ones onto this color so we don't have a gap
        if ( $states_handled == $state_count && $width_used < $progress_image_width ) {
            $width += $progress_image_width - $width_used;
        }

        push @state_elements, { state => $status, count => $states{$status}, width => $width };
    }

    if ($state_count > 1) {
        $has_multiple_states = 1;
    }
    
    ## build the messages line, if there are any
    my $messages_line = '';
    if (scalar @messages) {

        for my $msg (@messages) {
            ## show a counter if a message happened more than once
            if ($message_counts{$msg} > 1) {
                $messages_line .= "$msg <strong>($message_counts{$msg} times)</strong><br>";
            } else {
                $messages_line .= "$msg<br>";
            }
        }
    }
    
    ## "iterator" just means its on the distributed step.  make it more descriptive
    if ( $current_step =~ /^iterator /i ) {
        $current_step = "running distributed jobs ($current_step)";
    }
    
    ## we can adjust the default update interval here depending on what
    ##  state the component is in
    if ($component_state eq 'complete') {
        $update_interval = 0;
    } elsif ($component_state eq 'running') {
        $update_interval = 31;
    }

} else {
    my $dname = dirname("$pipeline");
    #Check for misconfigured pipeline
    if (-z $pipeline || -z "$pipeline.gz") {

	#Error creating XML file
	if(-s "$dname/replace_template_keys.stderr"){
	    $messages_line .= "Possible mismatch between component template and configuration file";
	    $messages_line .= "<br>".`cat $dname/replace_template_keys.stderr`;
	    $messages_line .= "<br>Check <a href='view_formatted_ini_source.cgi?pipeline_id=$pipelineid&file=$component_conf_varreplaced'>$component_conf_varreplaced</a>";
	}
	else{

	}
	$component_state = 'component configuration error';

	push @state_elements, { state => $component_state, count => 1, width => $progress_image_width };  
    }
    else{
	#Error creating config file
	if(-s "$dname/replace_config_keys.stderr"){
	    $messages_line .= "<br>Check <a href='view_formatted_ini_source.cgi?pipeline_id=$pipelineid&file=$component_conf_nonvarreplaced'>$component_conf_nonvarreplaced</a>";
	    $component_state = 'component configuration error';
	    push @state_elements, { state => $component_state, count => 1, width => $progress_image_width };  
	}
	else{
	    ## if the component.xml doesn't exist, that file-based subflow step in the pipeline.xml
	    #    hasn't started yet.  pipeline_summary.cgi should take care of those levels.
	    ## all we have here is the component state
	    $component_state = 'incomplete';
	    push @state_elements, { state => $component_state, count => 1, width => $progress_image_width };        
	}
    }
}

$tmpl->param( ACTION_COUNT => $command_count );
$tmpl->param( COMPONENT_CONFIG => -e $component_conf_varreplaced ? $component_conf_varreplaced : $component_conf_nonvarreplaced );
$tmpl->param( COMPONENT_STATE => $component_state );
$tmpl->param( CURRENT_STEP => $current_step );
$tmpl->param( HAS_MULTIPLE_STATES => $has_multiple_states );
$tmpl->param( MESSAGES_LINE => $messages_line );
$tmpl->param( PARENT_PIPELINE => $parent_pipeline );
$tmpl->param( PARENT_IS_RUNNING => $p_pipeline_state eq "running" ? 1 : 0 );
$tmpl->param( PARENT_PIPELINE_STATE => $p_pipeline_state );
$tmpl->param( PIPELINE => $pipeline );
$tmpl->param( PIPELINE_ID => $pipelineid);
$tmpl->param( PIPELINE_EXISTS => $pipeline_exists );
$tmpl->param( RUNTIME => $runtime );
$tmpl->param( STATE_ELEMENTS => \@state_elements );
$tmpl->param( UL_ID => $ul_id );
$tmpl->param( UPDATE_INTERVAL => $update_interval );


print $tmpl->output;

exit(0);

sub parseCommandSet {
    my ($commandSet, $fileparsed) = @_;
    
    ## get the component state if it hasn't been defined yet
    if ( $component_state eq 'unknown' ) {
        my $state = $commandSet->first_child('state') || 0;
        if ( $state ) {
            $component_state = $state->text;
        }
    }
    
    ## we need to get the status counts.  iterate through the commands and
    ##  get the states for each.  distributed subflow groups will count as
    ##  one command each here (mostly for parsing speed purposes)
    foreach my $command ( $commandSet->children('command') ) {

        ## all commands should have a state already, so I'll not check here.  this may
        #   need to be changed.
        my $state = 'unknown';
        if ( $command->has_child('state') ) {
            $state = $command->first_child('state')->text();
        }

        ## increase the count for this state
        $states{ $state }++;            

        ## if the state is running note it
        if ( $state eq 'running' ) {
            $current_step = $command->first_child('name')->text();
        }

        $command_count++;
    }

    if ( $commandSet->first_child('state') ) {
        $state = $commandSet->first_child('state')->text();
    }

    ## we don't want this to happen within gN.xml, only the parent component.xml
    if (! $got_time_info) {
        ($start_time, $end_time, $runtime) = &time_info( $commandSet );
        $got_time_info++;
    }

    ## all iterative components will have at least one commandSet to parse (file-based subflow iN.xml)
    #   these will be linked via fileName elements withinin commandSets
    if ( $fileparsed =~ /component.xml/) {
        for my $subflowCommandSet ( $commandSet->children('commandSet') ) {
            if ( $subflowCommandSet->first_child('state')->text() eq 'running' ) {
                $current_step = $subflowCommandSet->first_child('name')->text();
            }
        
            ## this command set should contain a fileName element
            my $fileName = $subflowCommandSet->first_child("fileName") || 0;
            if ($fileName) {
                if (-e $fileName->text || -e $fileName->text . '.gz') {
                    parseIteratorFile($fileName->text);
                }
            }
        }
    }
    
    my $text = $commandSet->sprint;
    while ( $text =~ m|<message>(.*?)</message>|g) {
        my $msg = $1;

        ## don't do anything with empty messages
        next unless length $msg > 0;

        ## don't include the "Command set with name: ? finished" messages
        next if ($msg =~ /Command set with name\:.+?finished/i);
        ## don't include "Job terminated." messages
        next if ($msg =~ /^Job terminated\.$/i);
        ## don't include "command finished" messages
        next if ($msg eq 'command finished');

        ## can we simplify this message?
        if ($msg =~ /SystemCommandProcessor line \d+\. (.+)/) {
            $msg = $1;
        }

        $msg =~ s/\</\&lt\;/g;
        $msg =~ s/\>/\&gt\;/g;

        ## add this to the messages array if it hasn't been seen yet.
        push @messages, $msg unless $message_counts{$msg}++;
    }
}


sub parseIteratorFile {
    ## this file should be the component's iN.xml
    my ($iN_xml) = @_;
    
    my $iN_xml_fh;
    if ($iN_xml =~ /\.gz/) {
        open($iN_xml_fh, "<:gzip", "$iN_xml") || die "can't read $iN_xml: $!"; 
    } elsif ( ! -e $iN_xml && -e "$iN_xml.gz" ) {
        open($iN_xml_fh, "<:gzip", "$iN_xml.gz") || die "can't read $iN_xml: $!"; 
    } else {
        open($iN_xml_fh, "<$iN_xml") || die "can't read $iN_xml: $!";       
    }
    
    ## all the commandSet elements in the iN.xml file will reference groups (gN.xml) files
    #   via fileName elements
    my $twig = XML::Twig->new(
                    twig_roots => {
                        'commandSet/fileName' => sub {
                                                     parseGroupFile( $_[1]->text );
                                                     $_[0]->purge();
                                                 },
                        'commandSet/commandSet/state' => sub {
                                                   $states{ $_[1]->text }++;
                                                   $command_count++;
                                              },
                    },
               );
    $twig->parse($iN_xml_fh);
}


sub parseGroupFile {
    ## this file should be a single gN.xml file
    #   it will have both commands and commandSets.  the commands should simply
    #   be counted by the commandSets reference gN.iter.xml files via fileName elements
    my ($gN_xml) = @_;
    
    my $gN_xml_fh;
    if ($gN_xml =~ /\.gz/) {
        open($gN_xml_fh, "<:gzip", "$gN_xml") || die "can't read $gN_xml: $!"; 
    } elsif ( ! -e $gN_xml && -e "$gN_xml.gz" ) {
        open($gN_xml_fh, "<:gzip", "$gN_xml.gz") || die "can't read $gN_xml: $!"; 
    } else {
        open($gN_xml_fh, "<$gN_xml") || die "can't read $gN_xml: $!";       
    }
    
    my $twig;
    
    if ( $max_command_resolution ) {
    
        $twig = XML::Twig->new(
                    twig_roots => {
                        'command/state' => sub {
                                               $states{ $_[1]->text }++;
                                               $command_count++;
                                               $_[0]->purge();
                                           },
                        'commandSet/fileName' => sub {
                                                     if ( -e $_[1]->text ) {
                                                        parseGroupIterFile( $_[1]->text );
                                                     }
                                                     
                                                     $_[0]->purge();
                                                 },
                    },
         );
    } else { 
        
        $twig = XML::Twig->new(
                    twig_roots => {
                        'command/state' => sub {
                                               $states{ $_[1]->text }++;
                                               $command_count++;
                                               $_[0]->purge();
                                           },
                    },
         );
    }
    $twig->parse($gN_xml_fh);
}

sub parseGroupIterFile {
    ## this file will be a single gN.iter.xml file
    #   it is made up of only commandSets, each referencing one of the iterated files.
    #   for command counts, we need to only look at the status elements of each.
    my ($gN_iter_xml) = @_;
    
    my $gN_iter_xml_fh;
    if ($gN_iter_xml =~ /\.gz/) {
        open($gN_iter_xml_fh, "<:gzip", "$gN_iter_xml") || die "can't read $gN_iter_xml: $!"; 
    } elsif ( ! -e $gN_iter_xml && -e "$gN_iter_xml.gz" ) {
        open($gN_iter_xml_fh, "<:gzip", "$gN_iter_xml.gz") || die "can't read $gN_iter_xml: $!"; 
    } else {
        open($gN_iter_xml_fh, "<$gN_iter_xml") || die "can't read $gN_iter_xml: $!";       
    }
    
    return;
    
    my $twig = XML::Twig->new(
                    twig_roots => {
                        'commandSet/status' => sub {
                                                    my ($t, $elt) = @_;
                                                    
                                                    for my $status ( $elt->children ) {
                                                        my $gi = $status->gi;
                                                        
                                                        next if ( $gi eq 'total' || $gi eq 'message' );
                                                        if ( $status->text ) {
                                                            $states{ $status->gi } += $status->text;
                                                            $command_count += $status->text;
                                                        }
                                                    }
                                                    
                                                    $t->purge();
                                               },
                    },
               );
    $twig->parse($gN_iter_xml_fh);
}











