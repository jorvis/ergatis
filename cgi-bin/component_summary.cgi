#!/usr/local/bin/perl -w

use strict;
use CGI;
use CGI::Carp qw(fatalsToBrowser);
use Date::Manip;
use Ergatis::ConfigFile;
use HTML::Template;
use Monitor;
use POSIX;
use XML::Twig;

my $q = new CGI;

print $q->header( -type => 'text/html' );

my $tmpl = HTML::Template->new( filename => 'templates/component_summary.tmpl',
                                die_on_bad_params => 1,
                              );

## will be like:
## /usr/local/annotation/TGA1/workflow/runtime/split_fasta/29134_test2/pipeline.xml
my $pipeline = $q->param("pipeline") || die "pass pipeline";

## will be like:
## split_fasta.test2
my $ul_id = $q->param("ul_id") || die "pass ul_id";

## will be like:
## /usr/local/annotation/AA1/workflow/runtime/pipeline/29671/pipeline.xml.instance
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
##  like: /usr/local/scratch/annotation/EHA1/workflow/runtime/iprscan/30835_default/pipeline.xml
my ($root_dir, $project, $component, $token, $pipelineid, $component_conf_varreplaced, $component_conf_nonvarreplaced);
if ($pipeline =~ m|(.+/(.+)/workflow/runtime/(.+?)/(\d+)_(.+?))/pipeline.xml|) {
    $component_conf_nonvarreplaced = "$1/component.conf.bld.ini";
    $component_conf_varreplaced = "$1/pipeline.config";
    ($project, $component, $pipelineid, $token) = ($2, $3, $4, $5);
} else {
    print "invalid instance file path format";
    exit(1);
}

if (-e $pipeline || -e "$pipeline.gz") {
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
    
    ## "Run subflow" just means its on the distributed step.  make it more descriptive
    if ( $current_step eq 'Run subflow' ) {
        $current_step = 'running distributed jobs';
    }
    
    ## we can adjust the default update interval here depending on what
    ##  state the component is in
    if ($component_state eq 'complete') {
        $update_interval = 0;
    } elsif ($component_state eq 'running') {
        $update_interval = 31;
    }

} else {
    ## all we have here is the component state
    $component_state = 'incomplete';
    push @state_elements, { state => $component_state, count => 1, width => $progress_image_width };        

    ## if the pipeline.xml doesn't exist, we need to check for errors during the component creation.
    ## we can only do this if the user passed the parent pipeline
    if ($parent_pipeline) {
    
        my $parent_pipeline_fh;
        if ($parent_pipeline =~ /\.gz/) {
            open($parent_pipeline_fh, "<:gzip", "$parent_pipeline") || die "can't read $parent_pipeline: $!"; 
        } elsif ( -e "$parent_pipeline.gz" ) {
            open($parent_pipeline_fh, "<:gzip", "$parent_pipeline.gz") || die "can't read $parent_pipeline: $!"; 
        } else {
            open($parent_pipeline_fh, "<$parent_pipeline") || die "can't read $parent_pipeline: $!";       
        }
    
        my $twig = new XML::Twig;
           $twig->parse($parent_pipeline_fh);
        my $commandSetRoot = $twig->root;
        my $parentCommandSet = $commandSetRoot->first_child('commandSet');
        
        ## find the commandSet for this component
        foreach my $commandSet ( $parentCommandSet->children('commandSet') ) {
            ## have we found it?
            if ( $commandSet->first_child('configMapId')->text() eq "component_$ul_id" ) {
                ## if it has a status, see if there are any messages;
                if ( $commandSet->first_child('status') ) {
                    
                    if ( $commandSet->first_child('status')->first_child('message') ) {
                        $messages_line .= $commandSet->first_child('status')->first_child('message')->text();
                    }
                }
                
                last;
            }
        }
    }
}

## if the component hasn't started yet, only the component.conf.bld.ini file will exist.  
$tmpl->param( ACTION_COUNT => $command_count );
$tmpl->param( COMPONENT_CONFIG => -e $component_conf_varreplaced ? $component_conf_varreplaced : $component_conf_nonvarreplaced );
$tmpl->param( COMPONENT_STATE => $component_state );
$tmpl->param( CURRENT_STEP => $current_step );
$tmpl->param( HAS_MULTIPLE_STATES => $has_multiple_states );
$tmpl->param( MESSAGES_LINE => $messages_line );
$tmpl->param( PARENT_PIPELINE => $parent_pipeline );
$tmpl->param( PIPELINE => $pipeline );
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

        if ( $command->first_child('state') ) {
           my $state = $command->first_child('state')->text();
        
            ## increase the count for this state
            $states{ $state }++;            
            
            ## if the state is running note it
            if ( $state eq 'running' ) {
                $current_step = $command->first_child('name')->text();
            }
            
        } else {
            ## state may not have been created yet (this should be rare)
            $states{unknown}++;
        }

        $command_count++;
    }

    if ( $commandSet->first_child('state') ) {
        $state = $commandSet->first_child('state')->text();
    }

    ## we don't want this to happen within groups.xml, only the parent pipeline.xml
    if (! $got_time_info) {
        ($start_time, $end_time, $runtime) = &time_info( $commandSet );
        $got_time_info++;
    }

    ## all iterative components will have a commandSet to parse (file-based subflow)
    my $subflowCommandSet = $commandSet->first_child("commandSet") || 0;
    if ($subflowCommandSet) {
        ## this command set should contain a fileName element
        my $fileName = $subflowCommandSet->first_child("fileName") || 0;
        if ($fileName) {
            if (-e $fileName->text || -e $fileName->text . '.gz') {
                parseComponentSubflow($fileName->text);
            }
        }
    }
    
    ## this is a terrible way to do this, as it doubles the memory required for
    ##  the twig.  it's functional, but needs to be replaced.
    ## don't look into the groups.xml here (yet).
    if ($fileparsed !~ /groups.xml/ ) {
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
}


sub parseComponentSubflow {
    ## currently this file should be the component's groups.xml
    my ($groupsXML) = @_;
    
    my $groupsXML_fh;
    if ($groupsXML =~ /\.gz/) {
        open($groupsXML_fh, "<:gzip", "$groupsXML") || die "can't read $groupsXML: $!"; 
    } elsif ( ! -e $groupsXML && -e "$groupsXML.gz" ) {
        open($groupsXML_fh, "<:gzip", "$groupsXML.gz") || die "can't read $groupsXML: $!"; 
    } else {
        open($groupsXML_fh, "<$groupsXML") || die "can't read $groupsXML: $!";       
    }
    
    my $twig = new XML::Twig;
       $twig->parse($groupsXML_fh);
    my $commandSetRoot = $twig->root;
    
    parseCommandSet( $commandSetRoot->first_child('commandSet'), $groupsXML );
}

