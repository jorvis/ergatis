package Ergatis::Utils;

use strict;
use warnings;
use XML::Twig;
use Ergatis::Monitor;

my %component_list;
my $order = 0;

sub process_root {
    my ( $t, $e ) = @_;
    if ( $e->first_child_text('name') eq 'start pipeline:' ) {
        foreach my $child ( $e->children('commandSet') ) {
            process_child( $child, $e, 'null' );

            my $file = $child->first_child_text('fileName');

            #open XML file (regular or gzip) of commandSet
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

sub process_child {
    my ( $e, $parent, $component ) = @_;
    my $name;
	my $state;
    my $count = 1;

# Used to indicate everything above this level will be 'parallel', 'serial', or 'start pipeline'
    my $top = 0;

    if ( $component eq 'null' )
    {    # want to establish the main pipeline components
        if (
            $parent->first_child_text('name') eq 'serial'
            ||    #components have one of these 3 parent names
            $parent->first_child_text('name') eq 'parallel'
            || $parent->first_child_text('name') eq 'remote-serial'
            || $parent->first_child_text('name') eq 'start pipeline:'
          )
        {
            $name = $e->first_child_text('name');

            # Set to top for major pipeline component or serial/parallel name
            $top = 1;
            if (   $name ne 'parallel'
                && $name ne 'serial'
                && $name ne 'remote-serial' )
            {
# component name is assigned... every child of component will pass cpu times to this key
            
			if ( $e->has_child('state') ) {
                $state = $e->first_child_text('state');
            }
                $component = $name;
                $component_list{$component}{'order'} = $order++;
				$component_list{$component}{'state'} = $state;
            }
        }
    }

#This mostly mirrors command from process_root... only dealing more with component commandSets
    foreach my $child ( $e->children('commandSet') ) {
        unless ( $child->first_child_text('name') ) {
            $name = $count++;
        } else {
            $name = $child->first_child_text('name');
        }
        process_child( $child, $e, $component );
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

	# Create twig XML to populate hash of component information
    my $twig = XML::Twig->new( 'twig_handlers' => { 'commandSet' => \&process_root } );
    $twig->parsefile($pipeline_xml);

	# Get progress rate information
	my ($total_complete, $total) = get_progress_rate(%component_list);
	my $failed_components = find_failed_components(\%component_list);
	
    my @failure_info = [ $failed_components, $stderr_files, $total_complete, $total ];
    return \@failure_info;
}

# Name: get_progress_rate
# Purpose: Parse through hash of component information and return the number of completed components and total components
# Args: hashref of component list information.  Can be created via XML::Twig from report_failure_info()
# Returns: Two variables
###	1) Integer of total components with a 'complete' status
### 2) Integer of total components

sub get_progress_rate {
	my $component_href = shift;

	# parse components array and count number of elements and 'complete' elements	
	my $total = scalar keys %$component_href;
	my @complete = grep { $component_href->{$_}->{'state'} eq 'complete'} keys %$component_href;
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
