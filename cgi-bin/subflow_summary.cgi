#!/usr/local/bin/perl -w

use strict;
use CGI;
use CGI::Carp qw(fatalsToBrowser);
use Date::Manip;
use POSIX;
use XML::Twig;

my $q = new CGI;

print $q->header( -type => 'text/html' );

my $xml_input = $q->param("xml_input") || die "pass xml_input";
my $subflow_state = 'unknown';
my $subflow_start = '';
my $subflow_end   = '';

print "    <h1>analysis steps</h1>\n";

## make sure the file exists, else print message
if (! -e $xml_input) {
    print "<div class='command'>steps not yet available</div>\n";
    exit;
}

my $twig = XML::Twig->new( twig_roots => {
                                'commandSet/state'     => sub {
                                                            my ($t, $elt) = @_;
                                                            $subflow_state = $elt->text;
                                                          },
                                'commandSet/startTime' => sub {
                                                            my ($t, $elt) = @_;
                                                            $subflow_start = $elt->text;
                                                          },
                                'commandSet/endTime'   => sub {
                                                            my ($t, $elt) = @_;
                                                            $subflow_end = $elt->text;
                                                          },
                                'command'              => \&process_command,
                           },
                      );
$twig->parsefile($xml_input);



sub process_command {
    my ($twig, $command) = @_;
    
    my $name  = $command->first_child('name')->text;
    my $state = $command->first_child('state')->text;

    my ($start_time, $end_time, $run_time) = time_info($command);
    
    my $time_label;
    ## if still running, display the start time
    if ($state eq 'running') {
        $time_label = "start: $start_time";
        
    ## else display the elapsed time
    } else {
        $time_label = "time: $run_time";
    }

    print <<CommAnD;
    <div class='command'>
        <div class='leftside'>
            <img class='status' src='/cram/status_$state.png' title='$state' alt='$state'>
            $name
        </div>
        <div class='rightside'>
            <span class='minor'>$time_label</span>
            <span class='minor'>n/a</span>
        </div>
    </div>
CommAnD

    my $ret_value = 'unknown';
    my $message = '';
    
    ## if there is a status and a message, grab it
    if ( $command->first_child('status') ) {
        if ( $command->first_child('status')->first_child('retValue') ) {
            $ret_value = $command->first_child('status')->first_child('retValue')->text;
        }
        
        if ( $command->first_child('status')->first_child('message') ) {
            $message = $command->first_child('status')->first_child('message')->text;
        }
        
        if ( $message ) {
            print <<messageBLOCK;
        <div class='messageblock'>
            return value: $ret_value<br>
            message: $message
        </div>
messageBLOCK
        }
    }

    $twig->purge;
}


sub time_info {
    my $command = shift;
    
    ## make sure we can at least get start time
    if (! $command->first_child('startTime') ) {
        return ('unavailable', 'unavailable', 'unavailable');
    }
    
    my $state = $command->first_child('state')->text;
    my $start_time_obj = ParseDate( $command->first_child('startTime')->text );
    my $start_time     = UnixDate( $start_time_obj, "%c" );
    
    my ($end_time_obj, $end_time);
    ## end time may not exist (if running, for example)
    if ( $command->first_child('endTime') ) {
        $end_time_obj   = ParseDate( $command->first_child('endTime')->text );
        $end_time       = UnixDate( $end_time_obj, "%c" );
    }

    ## we can calculate runtime only if start and end time are known, or if start is known and state is running
    my $runtime = '?';
    if ($start_time_obj) {
        if ($end_time_obj) {
            $runtime = strftime( "%H hr %M min %S sec", reverse split(/:/, DateCalc($start_time_obj, $end_time_obj)) );
        } elsif ($state eq 'running') {
            $runtime = strftime( "%H hr %M min %S sec", reverse split(/:/, DateCalc("now", $start_time_obj)) ) . ' ...';
        }
    }
    
    ## if hours or minutes are 00, take them off
    $runtime =~ s/00 .+? //g;
    
    ## 00 seconds isn't possible
    $runtime = '&lt; 01 sec' if $runtime eq '00 sec';
    
    return ($start_time, $end_time, $runtime);
}








