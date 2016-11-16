use Date::Manip;
use File::stat;
use POSIX;
use Ergatis::ConfigFile;
use Ergatis::Common;

use strict;

sub component_count_hash {
    my $pipeline_file = shift;
    my %components;    # = ( 'wu-blastp' => {count => 5} );

    my $ifh;
    if ( $pipeline_file =~ /\.gz/ ) {
        open( $ifh, "<:gzip", "$pipeline_file" )
          || die "can't read $pipeline_file: $!";
    } else {
        open( $ifh, "<$pipeline_file" ) || die "can't read $pipeline_file: $!";
    }

    my $t = XML::Twig->new(
        twig_roots => {
            'commandSet' => sub {
                my ( $t, $elt ) = @_;

                if (   $elt->first_child('name')
                    && $elt->first_child('name')->text() =~ /^(.+?)\./ )
                {
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

sub component_info_aref {
    my $pipeline_file = shift;
    my @components;

    my $ifh;
    if ( $pipeline_file =~ /\.gz/ ) {
        open( $ifh, "<:gzip", "$pipeline_file" )
          || die "can't read $pipeline_file: $!";
    } else {
        open( $ifh, "<$pipeline_file" ) || die "can't read $pipeline_file: $!";
    }

    my $t = XML::Twig->new(
        twig_roots => {
            'commandSet' => sub {
                my ( $t, $elt ) = @_;

                if (   $elt->first_child('name')
                    && $elt->first_child('name')->text() =~ /^(.+?)\.(.+)/ )
                {
                    push @components, { name => $1, token => $2 };

                    if ( $elt->has_child('state') ) {
                        $components[-1]{state} =
                          $elt->first_child('state')->text;

                        if (   $components[-1]{state} eq 'error'
                            || $components[-1]{state} eq 'failed' )
                        {
                            $components[-1]{error_count}++;
                        }
                    }

                    (
                        $components[-1]{start_time},
                        $components[-1]{end_time},
                        $components[-1]{run_time}
                    ) = &time_info($elt);
                }
            },
        },
    );
    $t->parse($ifh);

    return \@components;
}

sub process_command {
    my ( $twig, $command ) = @_;

    my %cmd_props = (
        command_string => 'unknown',
        is_command     => 1,
        message        => '',
        name           => 'unknown',
        return_value   => 'unknown',
        stderr         => 'not defined',
        stdout         => 'not defined',
        id             => 'unknown',
        state          => 'unknown',
    );

    if ( $command->has_child('name') ) {
        $cmd_props{name} = $command->first_child('name')->text;
    }

    if ( $command->has_child('state') ) {
        $cmd_props{state} = $command->first_child('state')->text;
    }

    if ( $command->has_child('id') ) {
        $cmd_props{id} = $command->first_child('id')->text;
    }

    my $type = $command->first_child('type')->text;

    ( $cmd_props{start_time}, $cmd_props{end_time}, $cmd_props{run_time} ) =
      time_info($command);

    ## can we get a return value?
    if (   $command->first_child('status')
        && $command->first_child('status')->first_child('retValue') )
    {
        $cmd_props{return_value} =
          $command->first_child('status')->first_child('retValue')->text;
    }

    ## if there is a status and a message, grab it
    if (   $command->first_child('status')
        && $command->first_child('status')->first_child('message') )
    {
        $cmd_props{message} =
          $command->first_child('status')->first_child('message')->text;
    }

    ## can we build a command line string?

    ## the command itself is within an <executable> element
    if ( $command->has_child('executable') ) {
        $cmd_props{command_string} = $command->first_child('executable')->text;
    }

    my $command_args = '';
    my $arg          = '';
    for my $param ( $command->children('param') ) {
        my $key   = $param->first_child('key')->text;
        my $value = $param->first_child('value')->text;

        ## catch stderr and stdout
        if ( $key eq 'stdout' ) {
            $cmd_props{stdout} = $value;

            ## is this stderr?
        } elsif ( $key eq 'stderr' ) {
            $cmd_props{stderr} = $value;

            ## else it must be a parameter of the command
        } else {
            ## if the command type is RunUnixCommand and the key doesn't start with
            ##  the string '-', we need to add it.  this should be fixed later.
            ##  since workflow does it, we have to do it
            if ( $type eq 'RunUnixCommand' && $key !~ /^\-/ ) {
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
    $cmd_props{command_string} =
      "$cmd_props{command_string} $command_args $arg";

    $twig->purge;

    return %cmd_props;
}

sub quota_string {
    my $repository_root = shift;

    my $string = '';

    ## need to see if the user has this enabled
    my $ergatis_cfg = new Ergatis::ConfigFile( -file => "ergatis.ini" );

    if ( $ergatis_cfg->val( 'display_settings', 'enable_quota_lookup' ) ) {

        if ( $repository_root =~ m|^/usr/local/annotation/| ) {
            $string = `/usr/local/common/getquota -N $repository_root`;
            if ( $string =~ /(\d+)\s+(\d+)/ ) {
                my ( $limit, $used ) = ( $1, $2 );
                $string = sprintf( "%.1f", ( $used / $limit ) * 100 )
                  . "\% ($used KB of $limit KB used)";
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
    if ( !$command->first_child('startTime') ) {
        return ( 'unavailable', 'unavailable', 'unavailable' );
    }

    my $state          = $command->first_child('state')->text;
    my $start_time_obj = ParseDate( $command->first_child('startTime')->text );
    my $start_time     = UnixDate( $start_time_obj, "%c" );

    my ( $end_time_obj, $end_time );
    ## end time may not exist (if running, for example)
    if ( $command->first_child('endTime') ) {
        $end_time_obj = ParseDate( $command->first_child('endTime')->text );
        $end_time = UnixDate( $end_time_obj, "%c" );
    }

    ## we can calculate runtime only if start and end time are known, or if start is known and state is running
    my $runtime = '?';

    ## doing it here manually because strftime was behaving badly (or I coded it badly)
    if ($start_time_obj) {
        my $diffstring;
        $runtime = '';

        if ($end_time_obj) {
            $diffstring = DateCalc( $start_time_obj, $end_time_obj );
        } else {
            $diffstring = DateCalc( $start_time_obj, "now" );
        }

        ## take out any non \d: characters
        $diffstring =~ s/[^0-9\:]//g;

        my @parts = split( /:/, $diffstring );

        ## years + months + weeks + days
        my $val =
          ( $parts[0] * 365 ) +
          ( $parts[1] * 30 ) +
          ( $parts[2] * 7 ) +
          ( $parts[3] );
        if ( $val > 1 ) {
            $runtime .= "$val days ";
        } elsif ( $val == 1 ) {
            $runtime .= "$val day ";
        }

        $runtime .= "$parts[4] hr "  if $parts[4];
        $runtime .= "$parts[5] min " if $parts[5];
        $runtime .= "$parts[6] sec";
    }

    $runtime = '&lt; 1 sec' if $runtime eq '0 sec';

    return ( $start_time, $end_time, $runtime );
}

###
# Retrieves all pipeline XML for the given project and some preliminary metadata
# on each pipeline to aid in pagination.
###
sub get_pipeline_quickstats {
    my ( $ergatis_cfg, $rdh, $pipeline_root ) = @_;
    my %pipeline_quickstats;

    my $per_account_pipelines =
      $ergatis_cfg->val( 'authentication', 'per_account_pipeline_security' );
    my $account_pipelines = get_account_pipelines($ergatis_cfg);

    foreach my $pipeline_id ( readdir $rdh ) {
        next unless ( $pipeline_id =~ /^[A-Z0-9]+$/ );
        next
          if ( $per_account_pipelines
            && !is_admin_user($ergatis_cfg)
            && !exists( $account_pipelines->{$pipeline_id} ) );

        my $pipeline_file =
          "$pipeline_root/$pipeline_id/pipeline.xml";   ## may be modified below

        ## if only the pipeline.xml exists, we can do less
        if ( !-e $pipeline_file ) {

            if ( -e "$pipeline_file.gz" ) {
                $pipeline_file .= '.gz';
            } else {
                ## didn't actually find a pipeline file
                next;
            }
        }

        my $filestat      = stat($pipeline_file);
        my $pipeline_user = getpwuid( $filestat->uid );
        my $last_mod      = $filestat->mtime;

        $pipeline_quickstats{$pipeline_id} = {
            pipeline_id => $pipeline_id
            ,    ## set again so we can directly export to HTML::Template
            last_mod      => $last_mod,
            path          => $pipeline_file,
            pipeline_user => $pipeline_user,
        };
    }

    return %pipeline_quickstats;
}

###
# The time_info method implemented using the XML::LibXML module
###
sub time_info_libxml {
    my $command = shift;

    my $start_time_str = $command->findvalue('startTime');
    if ( !$start_time_str ) {
        return ( 'unavailable', 'unavailable', 'unavailable' );
    }

    my $start_time_obj = ParseDate($start_time_str);
    my $start_time = UnixDate( $start_time_obj, "%c" );

    my ( $end_time_obj, $end_time );
    my $end_time_str = $command->findvalue('endTime');
    if ($end_time_str) {
        $end_time_obj = ParseDate($end_time_str);
        $end_time = UnixDate( $end_time_obj, "%c" );
    }

    my $runtime = "?";

    if ($start_time_obj) {
        my $diffstring;
        $runtime = '';

        if ($end_time_obj) {
            $diffstring = DateCalc( $start_time_obj, $end_time_obj );
        } else {
            $diffstring = DateCalc( $start_time_obj, "now" );
        }

        ## take out any non \d: characters
        $diffstring =~ s/[^0-9\:]//g;

        my @parts = split( /:/, $diffstring );

        ## years + months + weeks + days
        my $val =
          ( $parts[0] * 365 ) +
          ( $parts[1] * 30 ) +
          ( $parts[2] * 7 ) +
          ( $parts[3] );
        if ( $val > 1 ) {
            $runtime .= "$val days ";
        } elsif ( $val == 1 ) {
            $runtime .= "$val day ";
        }

        $runtime .= "$parts[4] hr "  if $parts[4];
        $runtime .= "$parts[5] min " if $parts[5];
        $runtime .= "$parts[6] sec";
    }

    $runtime = '&lt; 1 sec' if $runtime eq '0 sec';

    return ( $start_time, $end_time, $runtime );
}

###
# The component_count_aref subroutine implemented using XML::LibXML
##
sub get_component_hash {
    my $command = shift;
    my @components;

    foreach my $node ( $command->findnodes('//commandSet') ) {
        my $name = $node->findvalue('name');
        if ( $name && $name =~ /^(.+?)\.(.+)/ ) {
            push( @components, { name => $1, token => $2, error_count => 0 } );

            my $state = $node->findvalue('state');
            if ($state) {
                $components[-1]{'state'} = $state;

                if ( $state eq 'error' || $state eq 'failed' ) {
                    $components[-1]{error_count}++;
                }
            }

            (
                $components[-1]{'start_time'},
                $components[-1]{'end_time'},
                $components[-1]{'run_time'}
            ) = time_info_libxml($node);
        }
    }

    return \@components;
}

###
# The component_count_hash subroutine implemented in XML::LibXML
###
sub get_component_list {
    my $command = shift;
    my %components;
    my $component_aref  = [];
    my $component_count = 0;

    foreach my $node ( $command->findnodes('//commandSet') ) {
        my $name = $node->findvalue('name');
        if ( $name && $name =~ /^(.+?)\./ ) {
            $components{$1}{'count'}++;

            my $state = $node->findvalue('state');
            if ($state) {
                if ( $state eq 'error' || $state eq 'failed' ) {
                    $components{$1}{'error_count'}++;
                }
            }
        }
    }

    # HTML::Template requires us to pass any complex data structures in an
    # array of hashes so we need to convert
    foreach my $component ( sort keys %components ) {
        $component_count += $components{$component}{'count'};

        push(
            @{$component_aref},
            {
                'name'        => $component,
                'count'       => $components{$component}{'count'},
                'error_count' => $components{$component}{'error_count'} || 0
            }
        );
    }

    return ( $component_count, $component_aref );
}

###
# Generates the range of pipeline ID's that should be displayed for this
# specific page of pipelines.
###
sub get_pipeline_ids_range {
    my ( $view, $min_pipeline_pos, $max_pipeline_pos, $pipelines_per_page,
        $quickstats_ref )
      = @_;
    my %pipeline_quickstats = %$quickstats_ref;
    my @pipelines_to_parse  = ();

    if ( $view eq 'group' || $view eq 'component' ) {
        ## No current pagination on grouped or component-based pipeline list views
        @pipelines_to_parse = keys %pipeline_quickstats;
    } else {
        ## Sort and record which to parse
        my $current_pipeline_pos = 0;
        for my $pipeline_id (
            sort {
                $pipeline_quickstats{$b}{last_mod}
                  cmp $pipeline_quickstats{$a}{last_mod}
            } keys %pipeline_quickstats
          )
        {
            $current_pipeline_pos++;

            ## Check and see if this pipeline is in our page range
            next
              if ( $current_pipeline_pos < $min_pipeline_pos
                || $current_pipeline_pos > $max_pipeline_pos );

            push( @pipelines_to_parse, $pipeline_id );

            if ( scalar(@pipelines_to_parse) >= $pipelines_per_page ) {
                last;
            }
        }
    }

    return @pipelines_to_parse;
}

###
# Calculates several statistics/values that are needed to properly pagination
# the pipelines for this project.
###
sub prepare_pipeline_pagination {
    my (
        $view_page_num,      $repository_root, $max_pages_to_show,
        $pipelines_per_page, $pipeline_count
    ) = @_;
    my $page_count;

    if ( $pipeline_count > $pipelines_per_page ) {
        if ( $pipeline_count % $pipelines_per_page ) {
            $page_count = int( $pipeline_count / $pipelines_per_page ) + 1;
        } else {
            $page_count = int( $pipeline_count / $pipelines_per_page );
        }
    } else {
        $page_count = 1;
    }

    ## Calculate the min/max pipeline positions to show based on the page requested
    my $min_pipeline_pos;
    my $max_pipeline_pos;

    if ( $view_page_num == 1 ) {
        $min_pipeline_pos = 1;
    } else {
        $min_pipeline_pos =
          ( $pipelines_per_page * ( $view_page_num - 1 ) ) + 1;
    }

    if ( $view_page_num * $pipelines_per_page >= $pipeline_count ) {
        $max_pipeline_pos = $pipeline_count;
    } else {
        $max_pipeline_pos = $view_page_num * $pipelines_per_page;
    }

    my $page_links = [];

    ## We can't show ALL pages since there could be dozens.
    ## Calculate the next highest/lowest interval of '10', then hide anything outside of that range
    my ( $next_lower_interval, $next_higher_interval );

    if ( $view_page_num == $max_pages_to_show ) {
        $next_lower_interval = ( $view_page_num / $max_pages_to_show );
    } elsif ( $view_page_num % $max_pages_to_show ) {
        $next_lower_interval =
          floor( $view_page_num / $max_pages_to_show ) * $max_pages_to_show + 1;
    } else {
        $next_lower_interval = $view_page_num - $max_pages_to_show + 1;
    }

    $next_higher_interval =
      ceil( $view_page_num / $max_pages_to_show ) * $max_pages_to_show + 1;

    ## These are calculated later and control the display of the boxes with '...' on either
    ## of the pagination ranges
    my $show_pre_continuation  = $next_lower_interval > 1            ? 1 : 0;
    my $show_post_continuation = $next_higher_interval < $page_count ? 1 : 0;

    for ( my $i = 1; $i <= $page_count; $i++ ) {
        next if $i >= $next_higher_interval || $i < $next_lower_interval;

        push @$page_links,
          {
            page_num  => $i,
            is_active => $i == $view_page_num ? 1 : 0,
            url =>
              "./pipeline_list.cgi?repository_root=$repository_root&page_num=$i"
          };
    }

    return ( $page_count, $min_pipeline_pos, $max_pipeline_pos,
        $show_pre_continuation, $show_post_continuation, $page_links );
}

1 == 1;
