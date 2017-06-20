package Ergatis::Utils;

use strict;
use warnings;
use Date::Format;
use Date::Parse;
use XML::Twig;
use Module::Load;

use Exporter qw(import);

our @EXPORT_OK = qw(build_twig create_progress_bar update_progress_bar handle_component_status_changes report_failure_info);

### GLOBAL
my %component_list;
my $order;


# Name: build_twig
# Purpose: Start building the XML Twig for the pipeline
# Args: path to a pipeline XML file
# Returns: Hashref of the components and properties

sub build_twig {
    my $pipeline_xml = shift;
    # Reset the order variable and component list hashes
    $order = 0;
    %component_list = ();	# Will be filled in at process_child subroutine
    # Create twig XML to populate hash of component information
    my $twig = XML::Twig->new( 'twig_handlers' => { 'commandSet' => \&process_root } );
    $twig->parsefile($pipeline_xml);
    return \%component_list;
}

# Name: process_root
# Purpose: Process the base pipeline XML using XML:Twig.
# Args: XML Twig and Twig element objects.  Both are passed as args, when creating a subroutine for the 'twig_handlers' property of a new XML::Twig object
# Returns: Nothing

sub process_root {
    my ( $t, $e ) = @_;
    # The highest commandSet tag is the 'start pipeline' tag
    if ( $e->first_child_text('name') =~ /start pipeline/ ) {
        # This will process all XML under the element's 'commandSet' child tag.
        foreach my $child ( $e->children('commandSet') ) {
            process_child( $child, $e, 'null' );

            #open XML file (regular or gzip) of commandSet if 'filename' tag exists
            my $file = $child->first_child_text('fileName');
            my $fh = open_fh($file);
            if ( defined $fh ) {
                my $twig = XML::Twig->new(
                    'twig_handlers' => {
                        'commandSetRoot' => sub {
                            my ( $t, $e ) = @_;
                            foreach my $ch ( $e->children('commandSet') ) {
                                process_child( $ch, $e, 'null' );
                            }
                        }
                    }
                );
                $twig->parse($fh);
            }
        }
    }
}

# Name: process_child
# Purpose: Process a nested pipeline XML using XML:Twig.  Keeps track of the order of appearance for each component as well as the running status of that component
# Args: Three arguments
### 1 - Child <commandSet> element in an XML
### 2 - Parent element to child element from first argument
### 3 - Component name associated with parent element.  If arg is 'null', then it is assumed that the parent XML element was either a group iteration number or a Workflow command (i.e. serial or parallel)
# Returns: Nothing

sub process_child {
    my ( $e, $parent, $component ) = @_;
    my $name;
    my $state;
    my $count = 1;

    # Handle cases where the XML passed is not for a particular component
    if ( $component eq 'null' ) {
        if (
            $parent->first_child_text('name') eq 'serial'
            ||    #components have one of these 3 parent names
            $parent->first_child_text('name') eq 'parallel'
            || $parent->first_child_text('name') eq 'remote-serial'
            || $parent->first_child_text('name') =~ /start pipeline/
          )
        {
            $name = $e->first_child_text('name');

            if (   $name ne 'parallel'
                && $name ne 'serial'
                && $name ne 'remote-serial' )
            {
# component name is assigned... every child of component will pass cpu times to this key

                if ( $e->has_child('state') ) {
                    $state = $e->first_child_text('state');
                }
                $component = $name;
                # Start and End times are in ISO-8601 format
                my $start   = $e->first_child_text('startTime') if $e->has_child('startTime');
                my $end     = $e->first_child_text('endTime') if $e->has_child('endTime');
                if (defined $start && defined $end) {
                    my $elapsed_str = get_elapsed_time( $start, $end );
                    $component_list{$component}{'wall'} = $elapsed_str;
                }
                $component_list{$component}{'order'} = $order++;
                $component_list{$component}{'state'} = $state;
                @{$component_list{$component}{'command'}} = ();
            }
        }
    } else {
        # Perform actions on identified components
		process_failed_commands($e, $component);
    }

#This mostly mirrors command from process_root... only dealing more with component commandSets
    foreach my $child ( $e->children('commandSet') ) {
        process_child( $child, $e, $component );

        # If there is an XML file, set the file handle for opening and process that
        my $file = $child->first_child_text('fileName');
        my $fh   = open_fh($file);
        if ( defined $fh ) {
            my $twig = XML::Twig->new(
                'twig_handlers' => {
                    'commandSetRoot' => sub {
                        my ( $t, $e ) = @_;
                        foreach my $ch ( $e->children('commandSet') ) {
                            process_child( $ch, $e, $component );
                        }
                    }
                }
            );
            $twig->parse($fh);
        }
    }
}

# Name: process_failed_commands
# Purpose: Iterate through the component's nested XML and add stderr files and command strings to the component hash if the command failed
# Args: XML Twig element, and the component to store stderr data for
# Returns: Nothing 
sub process_failed_commands {
    my ($e, $component) = @_;
    if ($component_list{$component}{'state'} =~ /(failed|error)/){
        foreach my $command_child ( $e->children('command') ) {
            if ($command_child->has_child('state') 
              && $command_child->first_child_text('state') =~ /(failed|error)/) {
			  	my $command_string = '';
				my $command_args = '';
				my $arg          = '';
				my $stdout = '';
				my $stderr = '';

				## the command itself is within an <executable> element
				if ( $command_child->has_child('executable') ) {
					$command_string = $command_child->first_child('executable')->text;
				} else { 
					# If there is no executable just skip parsing the rest of this command
					next; 
				}
                foreach my $param_child ( $command_child->children('param') ) {
					my $key   = $param_child->first_child('key')->text;
					my $value = $param_child->first_child('value')->text;

					## catch stderr and stdout
					if ( $key eq 'stdout' ) {
						$stdout = $value;
					} elsif ( $key eq 'stderr' ) {
						$stderr = $value;
					} else {
						## if the command type is RunUnixCommand and the key doesn't start with
						##  the string '-', we need to add it.  this should be fixed later.
						##  since workflow does it, we have to do it

					    my $type = $command_child->first_child('type')->text;
						if ( $type eq 'RunUnixCommand' && $key !~ /^\-/ ) {
							$key = '--' . $key;
						}
						$command_args .= " $key=$value";
					}
				}

				## snatch the arg element if there was one
				if ( $command_child->first_child('arg') ) {
					$arg = $command_child->first_child('arg')->text;
				}

				## finish the command string build
				$command_string = "$command_string $command_args $arg";
				$command_string .= " 1>$stdout" if (length $stdout);
				$command_string .= " 2>$stderr" if (length $stderr);

            	push @{$component_list{$component}{'command'}}, [$command_string, $stderr];
            }
        }
    }
}


# Name: report_failure_info
# Purpose: Gather information to help diagnose a failure in the pipeline.
# Args: Pipeline XML file, (optional) pipeline ID
# Returns: An hash reference with the following elements:
### Array reference of failed components
### ID of the pipeline
### Hashref of component keys, for each key an array of stderr files
### Hashref of component keys, for each key an array of command_strings
### Number of components that have completed running
### Total number of components in the pipeline

sub report_failure_info {
    my $pipeline_xml = shift;
	my $pipeline_id = shift;

	# If pipeline ID was not passed, define it here
	if (! defined($pipeline_id) && $pipeline_xml =~ /runtime\/pipeline\/(\d+)\/pipeline\.xml/) {
		$pipeline_id = $1;
	}

    build_twig($pipeline_xml);
    # Get progress rate information
    my ($total_complete, $total) = get_progress_rate_from_href(\%component_list);
    my $failed_components = find_failed_components(\%component_list);
	my $cmd_info_href = get_command_info(\%component_list);

    my %failure_info = ('components' => $failed_components, 
	  'pipeline_id' => $pipeline_id,
	  'command_info' => $cmd_info_href,
      'complete_components' => $total_complete, 
      'total_components' => $total );
    return \%failure_info;
}

# Name: get_progress_rate_from_href
# Purpose: Parse through hash of component information and return the number of completed components and total components
# Args: hashref of component list information.  Can be created via XML::Twig from report_failure_info()
# Returns: Two variables
### 1) Integer of total components with a 'complete' status
### 2) Integer of total components

sub get_progress_rate_from_href {
	my $component_href = shift;

	# parse components array and count number of elements and 'complete' elements
	my $total = scalar keys %$component_href;
	my @complete = grep { $component_href->{$_}->{'state'} eq 'complete'} keys %$component_href;
	my $complete = scalar @complete;
	return ($complete, $total);
}

# Name: get_progress_rate_from_xml
# Purpose: Parse pipeline XML and return the number of completed components and total components
# Args: path to a pipeline XML file
# Returns: Two variables
### 1) Integer of total components with a 'complete' status
### 2) Integer of total components

sub get_progress_rate_from_xml {
    my $pipeline_xml = shift;
    build_twig($pipeline_xml);
    return get_progress_rate_from_href(\%component_list);
}

# Name: find_failed_components
# Purpose: Use a hash of component information and return all keys whose state is "error" or "failed"
# Args: Hashref of component list information.  Can be created via XML::Twig from report_failure_info()
# Returns: Array reference of failed components
sub find_failed_components {
    my $component_href = shift;
    my @failed_components = grep { $component_href->{$_}->{'state'} =~ /(failed|error)/ } (sort { $component_href->{$a}->{'order'} <=> $component_href->{$b}->{'order'} } keys %$component_href);
    return \@failed_components;
}

# Name: get_command_info
# Purpose: Compile a list of command_information for every component command in failed or error state
# Args: Hashref of component list information.  Can be created via XML::Twig from report_failure_info()
# Returns: Hashref listing the component and within each component key, a 2-dimensional array of command_information for each command.
# 		   1st dimension is each failed command for the component
# 		   2nd dimension is the failed command string and the stderr file path
sub get_command_info {
    my $component_href = shift;
    my %cmd_str;
    foreach my $component (sort { $component_href->{$a}->{'order'} <=> $component_href->{$b}->{'order'} } keys %$component_href) {
        $cmd_str{$component} = $component_href->{$component}->{'command'} if (scalar @{$component_href->{$component}->{'command'}} > 0);
    }
    return \%cmd_str;
}
# Name: handle_component_status_changes
# Purpose: Determine if any pipeline components have changed status since the last cycle, and handle accordingly
### Print to STDOUT if the component has entered 'running' status
### Print to STDOUT the elapsed time if a component has entered 'complete' status (may add failed statuses later)
# Args: Hashref of component data, and a hashref of component data before the latest XML:Twig build
# Returns: An arrayref of updated running components
sub handle_component_status_changes {
    my ($component_href, $prev_states) = @_;

    foreach my $component (
        sort { $component_href->{$a}->{'order'} <=> $component_href->{$b}->{'order'} } keys %$component_href
    ) {
        my $old_state = $prev_states->{$component};
        my $new_state = $component_href->{$component}->{'state'};
        next if ( $old_state eq "complete" ); # Complete components do not change
        next if ( $old_state eq $new_state ); # Skip unchanged components
        # Capitalize latest state of component
        my $printed;
        ($printed = $new_state) =~ s/([\w']+)/\u\L$1/g;
        # Handle the various updated component states

        # First a special case, where the component started and finished before the next sleep cycle
        if ($old_state eq "incomplete" && $new_state =~ /^(complete|error|failed)/) {
            print STDOUT "== Running: $component\n";
        }

        if ($new_state eq "running") {
            print STDOUT "== $printed: $component\n";
        } elsif ($new_state =~ /^(complete|error|failed)/) {
            my $elapsed = $component_href->{$component}->{'wall'};
            print STDOUT "==== $printed: $component $elapsed\n\n";
        }
    }
}

# Name: create_progress_bar
# Purpose: Create a progress bar to visually track the progress of the pipeline
# Args: Path to pipeline xml, and ID of pipeline
# Returns: A Term::ProgressBar object

sub create_progress_bar {
    my ($pipeline_xml, $id) = @_;
    my ($complete, $total) = get_progress_rate_from_xml($pipeline_xml);
    # Create the progress bar, but also remove from terminal when finished
	my $module = "Term::ProgressBar";
	load $module or croak("Cannot load module $module.  Does it need to be installed with CPAN first?");
    my $p_bar = Term::ProgressBar->new ({name => "Pipeline $id",
                                         count => $total,
                                         remove => 1});
    return $p_bar;
}

# Name: update_progress_bar
# Purpose: Update an existing progress bar based on number of complete components
# Args: A Term::ProgressBar object
# Returns: Nothing
sub update_progress_bar {
    my ($p_bar, $component_href) = @_;
    my ($complete, $total) = get_progress_rate_from_href($component_href);
    $p_bar->update($complete);
    $p_bar->message("$complete out of $total components have completed");
}

# Name: get_elapsed_time
# Purpose: Get difference of two wall time periods using Date::Parse
# Args: A start time value and end time value
# Returns: String of hh:mm:dd for elapsed time
sub get_elapsed_time {
    my ( $start, $end ) = @_;
    my ( $s, $e, $elapsed );
    $s       = str2time($start);
    $e       = str2time($end);
    $elapsed = $e - $s;
    my $elapsed_string = sec2string($elapsed);
    #print "ELAPSED TIME IS: " . $elapsed_string, "\n";
    return $elapsed_string;
}

# Name: sec2string
# Purpose: Converts seconds to strings.  Taken from http://www.perlmonks.org/?node_id=30392 ('best' formula)
# Args: Time values in seconds
# Returns: String of time in hh:mm:ss
sub sec2string {
    my $s = shift;

    return sprintf "00:00:%02d", $s if $s < 60;

    my $m = $s / 60; $s = $s % 60;
    return sprintf "00:%02d:%02d", $m, $s if $m < 60;

    my $h = $m /  60; $m %= 60;
    return sprintf "%02d:%02d:%02d", $h, $m, $s if $h < 24;

    # If the 'hours' slot overflows, your program is probably taking too long to run ;-)
}

# Name: open_fh
# Purpose: Open a file for reading.  If no file exists, just return an undefined file handle
# Args: Path to a compressed or uncompressed file
# Returns: The file handle
sub open_fh {
    my $file = shift;
    my $fh;

    return $fh if ( !-e $file );
    if ( $file =~ /\.gz$/ ) {
        open( $fh, "<:gzip", $file ) || die "Can't read file $file: $!";
    } else {
        open( $fh, "<$file" ) || die "Can't read file $file: $!";
    }
    return $fh;
}

1 == 1;
