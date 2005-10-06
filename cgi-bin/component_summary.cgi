#!/usr/local/bin/perl -w

use strict;
use CGI;
use CGI::Carp qw(fatalsToBrowser);
use XML::Twig;
use Date::Manip;
use Monitor;
use POSIX;

my $q = new CGI;

print $q->header( -type => 'text/html' );

## will be like:
## /usr/local/scratch/annotation/TGA1/Workflow/split_fasta/29134_test2/pipeline.xml
my $pipeline = $q->param("pipeline") || die "pass pipeline";

## will be like:
## split_fasta_test2
my $ul_id = $q->param("ul_id") || die "pass ul_id";

my $progress_image_width = 500;
my $component_state = 'unknown';
my $command_count = 0;
my %states;
my @messages;
my %message_counts;
#time variables
my ($start_time, $end_time, $lastmodtime, $state, $runtime) = ('n/a', 'n/a', '', 'unknown', 'n/a');
my $got_time_info = 0;

## give colors as rgb values or hexidecimal
my %colors = (
                complete    => 'rgb(0,200,0)',      ## green
                incomplete  => 'rgb(75,75,75)',     ## dark grey
                failed      => 'rgb(200,0,0)',      ## red
                pending     => 'rgb(200,200,0)',    ## yellow
                errors      => 'rgb(200,0,0)',      ## red
                error       => 'rgb(200,0,0)',      ## red
                running     => 'rgb(0,0,200)',      ## blue
                waiting     => 'rgb(200,200,0)',    ## yellow
                interrupted => 'rgb(200,0,200)',    ## purple
                total       => 'rgb(0,0,0)',        ## black 
             );

## we can parse some information out of the standardized instance file path
##  like: /usr/local/scratch/annotation/EHA1/Workflow/iprscan/30835_default/pipeline.xml
my ($root_dir, $project, $component, $token, $pipelineid, $component_conf);
if ($pipeline =~ m|(.+/(.+)/Workflow/(.+?)/(\d+)_(.+?))/pipeline.xml|) {
    #$component_conf = "$1/component.conf.bld.ini";
    $component_conf = "$1/pipeline.config";
    ($project, $component, $pipelineid, $token) = ($2, $3, $4, $5);
} else {
    print "invalid instance file path format";
    exit;
}

if (-e $pipeline) {

    my $twig = new XML::Twig;
       $twig->parsefile($pipeline);
    my $commandSetRoot = $twig->root;
    parseCommandSet( $commandSetRoot->first_child('commandSet'), $pipeline );

    ## build the line that lists each status and its count
    my $status_list_line = '';
    if (scalar keys %states > 1) {
        $status_list_line = '<li>states: ';
        
        for my $status (sort keys %states) {
            $status_list_line .= "$status (<span style='color:$colors{$status};'>$states{$status}</span>), ";
        }
        
        ## take off the trailing comma
        if ($status_list_line =~ /(.+)\,\s*$/) {
            $status_list_line = $1;
        }
        
        $status_list_line .= "</li>\n";
    }

    ## build the "image" div contents that represent the different states
    my $status_image = '';
    for my $status (sort keys %states) {
        ## each status gives a percentage of the total command_count
        $status_image .= '<div class="status_bar_portion" style="width: ' . 
                           int( ($states{$status} / $command_count) * $progress_image_width) . 'px; background-color: ' . 
                           ($colors{$status} || 'rgb(0,0,0)') . ';">' . "</div>\n";
    }

    ## build the messages line, if there are any
    my $messages_line = '';
    if (scalar @messages) {
        $messages_line = '<li class="messages"><b>messages:</b><br>';
        for my $msg (@messages) {
            ## show a counter if a message happened more than once
            if ($message_counts{$msg} > 1) {
                $messages_line .= "$msg <strong>($message_counts{$msg} times)</strong><br>";
            } else {
                $messages_line .= "$msg<br>";
            }
        }
        $messages_line .= "</li>\n";
    }

    ## we can adjust the default update interval here depending on what
    ##  state the component is in
    my $update_interval = 61;
    if ($component_state eq 'complete') {
        $update_interval = 0;
    } elsif ($component_state eq 'running') {
        $update_interval = 31;
    }

    ## print the component summary HTML
    print <<ComPONENTSummary;
    <h1><div class="component_label"><b>component</b>: $ul_id</div><div class="timer" id="${ul_id}_timer_label">update in <span id='${ul_id}_counter'>10</span>s</div></h1>
    <li><div class="component_progress_image">$status_image</div></li>
    <li>state: <span style='color: $colors{$component_state}'>$component_state</span> actions: $command_count</li>
    $status_list_line
    <b>runtime</b>: $runtime<br />
    $messages_line
    <li class="actions">
        <a href="./view_component.cgi?pipeline_xml=$pipeline"><img class='navbutton' src='/ergatis/button_blue_view.png' alt='view' title='view'></a>
        <a href="./view_formatted_xml_source.cgi?file=$pipeline" target="_blank"><img class='navbutton' src='/ergatis/button_blue_xml.png' alt='xml' title='xml'></a>
        <a href="./view_formatted_ini_source.cgi?file=$component_conf" target="_blank"><img class='navbutton' src='/ergatis/button_blue_config.png' alt='config' title='config'></a> 
        <a onclick="requestComponentUpdate('$pipeline', '$ul_id')"><img class='navbutton' src='/ergatis/button_blue_update.png' alt='update' title='update'></a>
        <a onclick="stopAutoUpdate('$ul_id')"><img class='navbutton' src='/ergatis/button_blue_stop_update.png' alt='stop update' title='stop update'></a>
    </li>
    <li class="pass_values">
        <span id="${ul_id}_continue_update">$update_interval</span>
    </li>
ComPONENTSummary

} else {
    ## print the incomplete component summary HTML
    printIncompleteSummary();
}

sub parseCommandSet {
    my ($commandSet, $fileparsed) = @_;
    
    ## get the component state if it hasn't been defined yet
    if ( $component_state eq 'unknown' ) {
        my $state = $commandSet->first_child('state') || 0;
        if ( $state ) {
            $component_state = $state->text;
        }
    }
    
    ## pull the status section for this commandset
    my $status = $commandSet->first_child('status') || 0;
    if ($status) {
        for my $type ($status->children) {
            ## don't do the 'total' or 'message' types here
            next if ($type->gi eq 'total' || $type->gi eq 'message');
            
            if ($type->text) {
                $states{$type->gi} += $type->text;
                $command_count += $type->text;
            }
        }
    }

    if ( $commandSet->first_child('state') ) {
        $state  = $commandSet->first_child('state')->text();
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
            if (-e $fileName->text) {
                parseComponentSubflow($fileName->text);
            }
        }
    }
    
    ## this is a terrible way to do this, as it doubles the memory required for
    ##  the twig.  it's functional, but needs to be replaced.
    ## don't look into the groups.xml here (yet).
    if ($fileparsed !~ /groups.xml$/) {
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
    
    my $twig = new XML::Twig;
       $twig->parsefile($groupsXML);
    my $commandSetRoot = $twig->root;
    parseCommandSet( $commandSetRoot->first_child('commandSet'), $groupsXML );

}




sub printIncompleteSummary {
    print <<ComponentNOTyetCreated;
    <h1><div class="component_label"><b>component</b>: $ul_id</div><div class="timer" id="${ul_id}_timer_label">update in <span id='${ul_id}_counter'>10</span>s</div></h1>
    <li><div class="component_progress_image"><div class="status_bar_portion" style="width: 500px; background-color: $colors{incomplete}"></div></div></li>
    <li>state: <span style='color: $colors{incomplete}'>incomplete</span></li>
    <li class="actions">
        <img class='navbutton' src='/ergatis/button_grey_view.png' alt='view' title='view'> 
        <img class='navbutton' src='/ergatis/button_grey_xml.png' alt='xml' title='xml'>  
        <img class='navbutton' src='/ergatis/button_grey_config.png' alt='config' title='config'>  
        <a onclick="requestComponentUpdate('$pipeline', '$ul_id')"><img class='navbutton' src='/ergatis/button_blue_update.png' alt='update' title='update'></a>
        <a onclick="stopAutoUpdate('$ul_id')"><img class='navbutton' src='/ergatis/button_blue_stop_update.png' alt='stop update' title='stop update'></a>
    </li>
    <li class="pass_values">
        <span id="${ul_id}_continue_update">60</span>
    </li>
ComponentNOTyetCreated
}

