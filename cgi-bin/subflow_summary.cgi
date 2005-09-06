#!/usr/local/bin/perl -w

use strict;
use CGI;
use CGI::Carp qw(fatalsToBrowser);
use Date::Manip;
use File::Basename;
use POSIX;
use XML::Twig;

my $q = new CGI;

print $q->header( -type => 'text/html' );

my $xml_input = $q->param("xml_input") || die "pass xml_input";
my $subflow_name  = basename($xml_input, '.xml');
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
                                                            print "<span class='hidden' id='${subflow_name}_state'>$subflow_state</span>\n";
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
    my $id    = $command->first_child('id')->text;
    my $type  = $command->first_child('type')->text;

    my ($start_time, $end_time, $run_time) = time_info($command);
    
    ## can we get a return value?
    my $return_value = 'unknown';
    if ( $command->first_child('status') && $command->first_child('status')->first_child('retValue') ) {
         $return_value = $command->first_child('status')->first_child('retValue')->text;
    }

    ## can we build a command line string?
    my $command_string = 'unknown';
    my $command_args = '';
    my $arg = '';
    my $stderr = 'not defined';
    my $stdout = 'not defined';
    for my $param ( $command->children('param') ) {
        my $key   = $param->first_child('key')->text;
        my $value = $param->first_child('value')->text;
        
        ## is this the command?
        if ($key eq 'command') {
            $command_string = $value;
        
        ## is this stdout?
        } elsif ($key eq 'stdout') {
            $stdout = $value;
            
        ## is this stderr?
        } elsif ($key eq 'stderr') {
            $stderr = $value;
            
        ## else it must be a parameter of the command
        } else {
            ## if the command type is RunUnixCommand and the key doesn't start with
            ##  the string '--', we need to add it.  this should be fixed later.
            ##  since workflow does it, we have to do it
            if ($type eq 'RunUnixCommand' && $key !~ /^\-\-/) {
                $key = '--' . $key;
            }
            
            $command_args .= " $key=$value";
        }
    }
    
    ## snatch the arg element if there was one
    if ( $command->first_child('arg') ) {
        $arg = $command->first_child('arg')->text;
    }
    
    ## finish the command string build
    $command_string = "$command_string $command_args $arg";

    print <<CommAnD;
    <div class='command'>
        <div class='leftside'>
            <img class='status' src='/cram/status_$state.png' title='$state' alt='$state'>
            $name
        </div>
        <div class='rightside'>
            <span class='minor'>$run_time</span>
            <span class='infolabel' id='${id}_infolabel'   onclick='toggle_cmd_info("$id")'>show info</span>
        </div>
    </div>
    <div class='cmdinfo' id='$id' style='display: none;'>
        <table>
            <tr>
                <th>workflow id:</th><td>$id</td>
            </tr>
            <tr>
                <th>state:</th><td>$state</td>
            </tr>
            <tr>
                <th>start time:</th><td>$start_time</td>
            </tr>
            <tr>
                <th>end time:</th><td>$end_time</td>
            </tr>
            <tr>
                <th>duration:</th><td>$run_time</td>
            </tr>
            <tr>
                <th>return value:</th><td>$return_value</td>
            </tr>
            <tr>
                <th>stdout:</th><td>$stdout</td>
            </tr>
            <tr>
                <th>stderr:</th><td>$stderr</td>
            </tr>
            <tr>
                <th colspan='2'>command:</th>
            </tr>
            <tr>
                <td colspan='2'>$command_string</td>
            </tr>
        </table>
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
    
    ## take off leading zero (if present)
    if ($runtime =~ /^0(.+)/) {
        $runtime = $1;
    }
    
    ## 00 seconds isn't possible
    $runtime = '&lt; 1 sec' if $runtime eq '0 sec';
    
    return ($start_time, $end_time, $runtime);
}








