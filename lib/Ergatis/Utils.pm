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
                $component = $name;
                $component_list{$component}{'Order'} = $order++;
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

sub report_failure_info {
    my ($pipeline_xml) = shift;

     #my $twig =
     #  XML::Twig->new( 'twig_handlers' => { 'commandSet' => \&process_root } );
     #$twig->parsefile($pipe_xml);

	my ($total_complete, $total) = get_progress_rate($pipeline_xml);

    my @failure_info = [ $failed_components, $stderr_msgs, $total_complete, $total ];
    return \@failure_info;
}

# Name: get_progress_rate
# Purpose: Parse through the pipeline XML and return the number of completed components and total components
# Args: Pipeline XML file
# Returns: Two variables
###	1) Integer of total components with a 'complete' status
### 2) Integer of total components

sub get_progress_rate {
	my $xml = shift;

    my $ifh;
    if ( $xml =~ /\.gz/ ) {
        open( $ifh, "<:gzip", "$xml" )
          || die "can't read $xml: $!";
    } else {
        open( $ifh, "<$xml" ) || die "can't read $xml: $!";
    }

	my @components;

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
                    }
                }
            },
        },
    );
    $t->parse($ifh);

	# parse components array and count number of elements and 'complete' elements	
	my $total = scalar @components;
	my @complete = grep { $_{state} eq 'complete'} @components;
	my $complete = scalar @complete;
	return ($complete, $total);

}

sub find_failed_components {

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
