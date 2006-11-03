use Date::Manip;
use File::stat;
use POSIX;
use Ergatis::ConfigFile;

use strict;


sub component_count_hash {
    my $pipeline_file = shift;
    my %components; # = ( 'wu-blastp' => {count => 5} );

    #return %components;

    my $ifh;
    if ($pipeline_file =~ /\.gz/) {
        open($ifh, "<:gzip", "$pipeline_file") || die "can't read $pipeline_file: $!"; 
    } else {
        open($ifh, "<$pipeline_file") || die "can't read $pipeline_file: $!";       
    }

    my $t = XML::Twig->new();
       $t->parse($ifh);

    my $first_cs = $t->root->first_child('commandSet') || die "failed to pull the first commandSet";
    
    ## each command set under the first is a component
    for my $cs ( $first_cs->children('commandSet') ) {
        if ( $cs->first_child('name')->text =~ /(.+?)\./ ) {
            my $component_name = $1;
            $components{$component_name}{count}++;
            
            if ( $cs->has_child('state') ) {
                my $state = $cs->first_child('state')->text;
                
                if ( $state eq 'error' || $state eq 'failed' ) {
                    $components{$component_name}{error_count}++;
                }
            }            
        }
    }
    
    return %components;
}


sub process_command {
    my ($twig, $command) = @_;
    
    my %cmd_props = ( 
                      command_string => 'unknown',
                      is_command => 1, 
                      return_value => 'unknown',
                      stderr => 'not defined',
                      stdout => 'not defined',
                      message => '',
                    );
    
    $cmd_props{name}  = $command->first_child('name')->text;
    $cmd_props{state} = $command->first_child('state')->text;
    $cmd_props{id}    = $command->first_child('id')->text;
    my $type  = $command->first_child('type')->text;

    ($cmd_props{start_time}, $cmd_props{end_time}, $cmd_props{run_time} ) = time_info($command);
    
    ## can we get a return value?
    if ( $command->first_child('status') && $command->first_child('status')->first_child('retValue') ) {
         $cmd_props{return_value} = $command->first_child('status')->first_child('retValue')->text;
    }

    ## if there is a status and a message, grab it
    if ( $command->first_child('status') && $command->first_child('status')->first_child('message') ) {
        $cmd_props{message} = $command->first_child('status')->first_child('message')->text;
    }

    ## can we build a command line string?
    my $command_args = '';
    my $arg = '';
    for my $param ( $command->children('param') ) {
        my $key   = $param->first_child('key')->text;
        my $value = $param->first_child('value')->text;
        
        ## is this the command?
        if ($key eq 'command') {
            $cmd_props{command_string} = $value;
        
        ## is this stdout?
        } elsif ($key eq 'stdout') {
            $cmd_props{stdout} = $value;
            
        ## is this stderr?
        } elsif ($key eq 'stderr') {
            $cmd_props{stderr} = $value;
            
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
    $cmd_props{command_string} = "$cmd_props{command_string} $command_args $arg";

    $twig->purge;
    
    return %cmd_props;
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
