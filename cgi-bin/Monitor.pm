use Date::Manip;
use File::stat;
use POSIX;
use Ergatis::ConfigFile;

use strict;


sub component_count_hash {
    my $pipeline_file = shift;
    my %components;

    my $ifh;
    if ($pipeline_file =~ /\.gz/) {
        open($ifh, "<:gzip", "$pipeline_file") || die "can't read $pipeline_file: $!"; 
    } else {
        open($ifh, "<$pipeline_file") || die "can't read $pipeline_file: $!";       
    }


    my $t = XML::Twig->new( twig_roots => {
                                'commandSet' => sub {
                                                      my ($t, $elt) = @_;
                                                      
                                                          if ($elt->first_child('configMapId')->text() =~ /^component_(.+?)\./) {
                                                              $components{$1}{count}++;
                                                              
                                                              if ( $elt->has_child('state') ) {
                                                                  my $state = $elt->first_child('state')->text;
                                                                  if ( $state eq 'error' || $state eq 'failed' ) {
                                                                      $components{$1}{error_count}++;
                                                                  }
                                                              }
                                                          }
                                                      },
                                          },
                          );
    $t->parse($ifh);
    
    return %components;
}


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
            ##  the string '-', we need to add it.  this should be fixed later.
            ##  since workflow does it, we have to do it
            if ($type eq 'RunUnixCommand' && $key !~ /^\-/) {
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
            <img class='status' src='/ergatis/status_$state.png' title='$state' alt='$state'>
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

sub quota_string {
    my $repository_root = shift;
    
    my $string = '';
    
    ## need to see if the user has this enabled
    my $ergatis_cfg = new Ergatis::ConfigFile( -file => "ergatis.ini" );
    
    if ( $ergatis_cfg->val('settings', 'enable_quota_lookup') ) {

        if ($repository_root =~ m|^/usr/local/annotation/|) {
            $string = `/usr/local/common/getquota -N $repository_root`;
            if ($string =~ /(\d+)\s+(\d+)/) {
                my ($limit, $used) = ($1, $2);
                $string = sprintf("%.1f", ($used/$limit) * 100) . "\% ($used KB of $limit KB used)";
            } else {
                $string = "error parsing quota information: $string";
            }
        } else {
            $string = 'unavailable (unknown project area)';
        }
    } else {
        $string = 'quota information currently disabled';
    }
    
    return $string;
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

    ## doing it here manually because strftime was behaving badly (or I coded it badly)
    if ($start_time_obj) {
        my $diffstring;
        $runtime = '';
        
        if ($end_time_obj) {
            $diffstring = DateCalc($start_time_obj, $end_time_obj, undef, 1);
        } else {
            $diffstring = DateCalc($start_time_obj, "now", undef, 1);
        }
        
        ## take out any non \d: characters
        $diffstring =~ s/[^0-9\:]//g;

        my @parts = split(/:/, $diffstring);
        
        ## years + months + weeks + days
        my $val = ($parts[0] * 365) + ($parts[1] * 30) + ($parts[2] * 7) + ($parts[3]);
        if ($val > 1) {
            $runtime .= "$val days ";
        } elsif ($val == 1) {
            $runtime .= "$val day ";
        }
        
        $runtime .= "$parts[4] hr " if $parts[4];
        $runtime .= "$parts[5] min " if $parts[5];
        $runtime .= "$parts[6] sec";
    }

    $runtime = '&lt; 1 sec' if $runtime eq '0 sec';

    return ($start_time, $end_time, $runtime);
}



1==1;
