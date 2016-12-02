package Ergatis::Utils;

use strict;
use warnings;
use Date::Format;
use Date::Parse;
use XML::Twig;
use Term::ProgressBar;

use Exporter qw(import);

our @EXPORT_OK = qw(build_twig create_progress_bar update_progress_bar handle_component_status_changes);

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
	%component_list = ();
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
    if ( $component eq 'null' )
    {
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
                    $component_list{$component}{'Wall'} = $elapsed_str;
                }
                $component_list{$component}{'order'} = $order++;
				$component_list{$component}{'state'} = $state;
            }
        }
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

# Name: report_failure_info
# Purpose: Gather information to help diagnose a failure in the pipeline.
# Args: Pipeline XML file
# Returns: An array reference with the following elements:
###	1) Array reference of failed components
### 2) Array reference of paths to stderr files
### 3) Number of components that have completed running
### 4) Total number of components in the pipeline

sub report_failure_info {
    my ($pipeline_xml) = shift;
    build_twig($pipeline_xml);
	# Get progress rate information
	my ($total_complete, $total) = get_progress_rate_from_href(\%component_list);
	my $failed_components = find_failed_components(\%component_list);

    my $stderr_files = undef;   # Placeholder for now

    my @failure_info = [ $failed_components, $stderr_files, $total_complete, $total ];
    return \@failure_info;
}

# Name: get_progress_rate_from_href
# Purpose: Parse through hash of component information and return the number of completed components and total components
# Args: hashref of component list information.  Can be created via XML::Twig from report_failure_info()
# Returns: Two variables
###	1) Integer of total components with a 'complete' status
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
###	1) Integer of total components with a 'complete' status
### 2) Integer of total components

sub get_progress_rate_from_xml {
    my $pipeline_xml = shift;
    build_twig($pipeline_xml);

	# parse components array and count number of elements and 'complete' elements
	my $total = scalar keys %component_list;
	my @complete = grep { $component_list{$_}->{'state'} eq 'complete'} keys %component_list;
	my $complete = scalar @complete;
	return ($complete, $total);
}

# Name: find_failed_components
# Purpose: Use a hash of component information and return all keys whose state is "error" or "failed"
# Args: Hashref of component list information.  Can be created via XML::Twig from report_failure_info()
# Returns: Array reference of failed components
sub find_failed_components {
	my $component_href = shift;
	my @failed_components = grep { $component_href->{$_}->{'state'} =~ /(failed|error)/ } keys %$component_href;
	return \@failed_components
}

sub get_failed_stderr {

}

# Name: handle_component_status_changes
# Purpose: Determine if any pipeline components have changed status since the last cycle, and handle accordingly
### Print to STDOUT if the component has entered 'running' status
### Print to STDOUT the elapsed time if a component has entered 'complete' status (may add failed statuses later)
# Args: Hashref of component data, and a hashref of component data before the latest XML:Twig build
# Returns: An arrayref of updated running components
sub handle_component_status_changes {
    my ($component_href, $prev_states) = @_;

    foreach my $component (keys %$prev_states) {
        my $old_state = $prev_states->{$component};
        my $new_state = $component_href->{$component}->{'state'};
        next if ( $old_state eq "complete" ); # Complete components do not change
        next if ( $old_state eq $new_state ); # Skip unchanged components
        # Capitalize latest state of component
        my $printed;
        ($printed = $new_state) =~ s/([\w']+)/\u\L$1/g;
        # Handle the various updated component states
        if ($new_state eq "running") {
            print STDOUT "== $printed: $component\n";
        } elsif ($new_state =~ /^(complete|error|failed)/) {
            my $elapsed = $component_href->{$component}->{'Wall'};
            print STDOUT "==== $printed: $elapsed\n\n";
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
