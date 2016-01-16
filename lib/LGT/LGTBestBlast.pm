
=head1 NAME

LGTBestBlast - Run blast and filter for best hits.

=head1 SYNOPSIS

Need to put something useful here

=head1 DESCRIPTION

A module to run BLAST or take existing blast -m8 output and filter it for
best hit. Also appends Lineage information onto each hit.

=head1 AUTHOR - David R. Riley

e-mail: driley@som.umaryland.edu

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

package LGT::LGTBestBlast;
use strict;
use warnings;
use File::Basename;
use Carp;
$Carp::MaxArgLen = 0;
$|               = 1;

# Globals
my $clones_by_clone_id = {};
my $trace_lookup       = {};
my $filter_hits        = [];
my $lineages           = [];
my $gi2tax;
my $overallfile;
my $out1file;
my $out2file;
my $out3;
my $trace_mapping_file;
my $FILTER_MIN_OVERLAP = 50;
my $BLAST_BIN = '/usr/local/bin/blastall';
my $BLAST_CMD_ARGS = ' -p blastn -e 1 -m8';
my $filter_lineage;
my $min_overlap;

=head2 &filterBlast

 Title   : filterBlast
 Usage   : my $bestBlast = LGTBestBlast::filterBlast({fasta => $fasta,...})
 Function: Returns the path to filtered blast reports
 Returns :  Hash ref. of filtered blast hits.

            list_file     => list of the files below (overall,out1,out2)
            overall_blast => Lin1 & Lin2 hits,
            out1file_file => lineage1 hits,
            out2file_file => lineage2 hits,

 Args    : A hash containing potentially several config options:

           fasta - Path to a fasta file to search
           blast - Path to an existing blast output
           db - One or more fasta files to use a host references
           output_dir - One or more fasta files to use as donor references
           blast_bin - path to blast (Can also contain some arguments)
           lineage1 - Lineage for donor
           lineage2 - Lineage for host
           gitaxon - A GiTaxon object

=cut

sub filterBlast {
    my ($class, $args) = @_;

# Default filter min overlap is 50 (minimum length to filter out overlapping reads
    $min_overlap =
      $args->{filter_min_overlap} ? $args->{filter_min_overlap} : $FILTER_MIN_OVERLAP;
    $filter_lineage = $args->{filter_lineage};

    # Figure out what inputs we have
    my $input;
    if ( $args->{fasta} ) {
        $input = $args->{fasta};
    }
    elsif ( $args->{blast} ) {
        $input = $args->{blast};
    }

    if ( !$input ) {
        confess "Need to provide a fasta or a blast file for blast\n";
    }

    # Initialize gi2taxon db.
    $gi2tax = $args->{gitaxon};
    if ( !$gi2tax ) { confess "Need to provide a gitaxon object\n"; }

    # Get the basename of the input file.
    my ( $name, $directories, $suffix ) = fileparse( $input, qr/\.[^.]*/ );

    # Open the overall file
    $overallfile = $args->{overalloutput};
    if ( !$args->{overalloutput} && $args->{output_dir} ) {
        $overallfile = $args->{output_dir} . "/$name\_overall.out";

        # print STDERR "Printing to $overallfile\n";
    }
    open $out3, ">$overallfile"
      or confess "Unable to open overall output $overallfile\n";

    my $trace_mapping_file = $args->{trace_mapping};

    &_read_map_file( $args->{trace_mapping} );
    &_init_lineages( $args, $name );

    $args->{blast_bin} = $BLAST_BIN if !defined ($args->{blast_bin});
    my $fh;
	# Use blast results if already provided, otherwise run blast to get results
    if ( $args->{blast} ) {
		print STDERR "Opening blast input.\n";
        open( $fh, "<$input" ) or confess "Unable to open $input\n";
    }
    elsif ( !$args->{blast} ) {
        open( $fh, "-|",
"$args->{blast_bin} $BLAST_CMD_ARGS -d $args->{db} -i $input"
          )
          or confess
          "Unable to run: $args->{blast_bin} on: $input with db: $args->{db}\n";
    }
    &_process_file($fh);
    my $list = &_create_list_of_outputs($args);

    foreach my $lineage (@$lineages) {
        close $lineage->{handle}
          ; ##  or confess "=== &LGTBestBlast - Can't close output filehandle because: $! ===\n";
    }

    return {
        list_file     => $list,
        overall_blast => $overallfile,
        out1file_file => $out1file,
        out2file_file => $out2file
    };
}

# Read in the trace mapping file if one was provided
sub _read_map_file {
    my $map_file = shift;
    if ($map_file) {
        print STDERR "Reading the map file\n";
        open MAP, "<$map_file" or confess "Unable to open $map_file\n";
        while (<MAP>) {
            my @fields = split;
            $trace_lookup->{ $fields[0] } = {
                'template_id' => $fields[2],
                'trace_end'   => $fields[1]
            };
            if ( !$clones_by_clone_id->{ $fields[2] } ) {
                $clones_by_clone_id->{ $fields[2] } = {
                    'forward' => {},
                    'reverse' => {}
                };
            }
        }
        close MAP;
        print STDERR "Done reading the map file\n";
    }
}

# Initialize the lineage hashes for the donor, the host, and overall
sub _init_lineages {
    my $args = shift;
    my $name = shift;
    if ( $args->{lineage1} && $args->{lineage2} ) {
        $out1file = $args->{output1};
        if ( !$args->{output1} && $args->{output_dir} ) {
            $out1file = $args->{output_dir} . "/$name\_lineage1.out";
        }
        open my $out1, ">$out1file"
          or confess "Unable to open lineage1 output $out1file\n";

        $out2file = $args->{output2};
        if ( !$args->{output1} && $args->{output_dir} ) {
            $out2file = $args->{output_dir} . "/$name\_lineage2.out";
        }
        open my $out2, ">$out2file"
          or confess "Unable to open lineage2 output $out2file\n";

        push( @$lineages,
            (
                {
                    'lineage'   => $args->{lineage1},
                    'tophit'    => 1,                  #$args->{lineage1tophit},
                    'best_e'    => 100,			# Dummy value so next will always be better
                    'id'        => '',
                    'handle'    => $out1,
                    'best_rows' => [],
                    'name'      => 'lineage1'
                },
                {
                    'lineage'   => $args->{lineage2},
                    'tophit'    => 1,                  #$args->{lineage2tophit},
                    'best_e'    => 100,
                    'id'        => '',
                    'handle'    => $out2,
                    'best_rows' => [],
                    'name'      => 'lineage2'
                }
            )
        );
    }
    push( @$lineages,
        {
            'lineage' => 'cellular organisms',    ## KBS 01.05.13
            'tophit'    => 1,
            'best_e'    => 100,
            'id'        => '',
            'handle'    => $out3,
            'best_rows' => [],
            'name'      => 'overall'
        }
    );
}

# Process the BLAST m8 output.
sub _process_file {
    my $fh = shift;
    use Data::Dumper;
    while (my $line = <$fh>) {
        chomp $line;
        my @new_fields = split( /\t/, $line );
        my $tax;
        my $found_tax = 0;

        # If we already have lineage info in here we'll not append it again
        if ( $new_fields[14] ) {
            $found_tax = 1;
            $tax       = {
                'taxon_id' => $new_fields[12],
                'lineage'  => $new_fields[14],
                'name'     => $new_fields[13]
            };
        }
        else {
            $tax = $gi2tax->getTaxon( $new_fields[1] );
        }

		# If taxon info was not already in m8 file, add it
		if (!$found_tax) {
			push (@new_fields, ($tax->{taxon_id}, $tax->{name}, $tax->{lineage} ) );
		}
		# Die if we are missing key info
		carp "Unable to find name or lineage for taxon_id $tax->{'taxon_id'} in trace $new_fields[1]\n"
			unless $tax->{name};
        carp "Unable to find taxon info for $new_fields[1]\n" unless $tax->{taxon_id};

        my $fields = &_add_trace_info( \@new_fields );

        my $done = 0;
        foreach my $lineage (@$lineages) {
            $done = &_process_line( $fields, $tax, $lineage );
			if ( $done == 1) {
				#&_append_hits( $lineage->{best_rows}, $lineage->{name} );
			}
        }

		# Reinitialize hash for the next m8 query ID
        if ($done) {
            $filter_hits = [];
        }
    }
    close $fh;

    # here we'll take care of the last trace in the file.
    foreach my $lineage (@$lineages) {
		&_append_hits( $lineage->{best_rows}, $lineage->{name} );
        &_print_hits($lineage);
    }
}

# Append finished query ID information into the clone and trace hash
### Shaun Adkins - 11/16/15 - I don't see where the clone hash or trace hash are used anywhere
sub _append_hits {
    my $hits      = shift;
    my $list_name = shift;

    my $template_id = $hits->[0]->[12];

    my $clone = $clones_by_clone_id->{$template_id};

    map {
        # Check what strand we're on
        if ( $_->[13] eq 'F' ) {
            # Create the trace entry
            if ( !$clone->{'forward'}->{ $_->[0] } ) {
                $clone->{'forward'}->{ $_->[0] } = {
                    'trace_id'   => $_->[0],
                    "$list_name" => []
                };
            }
            if ( !$clone->{'forward'}->{ $_->[0] }->{$list_name} ) {
                $clone->{'forward'}->{ $_->[0] }->{$list_name} = [];
            }
            push(
                @{ $clone->{'forward'}->{ $_->[0] }->{$list_name} },
                {
                    'accession' => $_->[1],
                    'pid'       => $_->[2],

                }
            );
        } elsif ( $_->[13] eq 'R' ) {

            # Create the trace entry
            if ( !$clone->{'reverse'}->{ $_->[0] } ) {
                $clone->{'reverse'}->{ $_->[0] } = {
                    'trace_id'   => $_->[0],
                    "$list_name" => []
                };
            }
            if ( !$clone->{'reverse'}->{ $_->[0] }->{$list_name} ) {
                $clone->{'reverse'}->{ $_->[0] }->{$list_name} = [];
            }
            push(
                @{ $clone->{'reverse'}->{ $_->[0] }->{$list_name} },
                {
                    'accession' => $_->[1],
                    'pid'       => $_->[2],

                }
            );
        }
    } @$hits;
}

sub _process_line {
    my ( $fields, $tax, $lineage ) = @_;

    my $finished_id = 0;

    # If we are still on the current query ID
	#print $lineage->{id} . " LINEAGE " . $lineage->{lineage}\n";
	#print $fields->[0] . " FIRST FIELD\n";
	
	# If our best lineage ID has been identified as the current query ID...
	if ( $lineage->{id} eq $fields->[0] ) {
        if ( $lineage->{tophit} ) {

            # Determining how to handle equal or better hits
            if ( $fields->[10] == $lineage->{best_e} ) {
				# Before updating lineage hash, make sure the hit lineage matches current lineage hash iteration
                if ( $tax->{lineage} =~ /$lineage->{lineage}/ ) {
					# If our hit is equally as good as our best, then allow for more than one best hit
                    push( @{ $lineage->{best_rows} }, $fields );
                    # print STDERR "Here with: $lineage->{lineage\n";
                }
            } elsif ( $fields->[10] < $lineage->{best_e} ) {
            	# If hit is lower than our best lineage row, it becomes the new best lineage
                if ( $tax->{lineage} =~ /$lineage->{lineage}/ ) {
                    $lineage->{best_e}    = $fields->[10];
                    $lineage->{best_rows} = [$fields];
                }
            }
        }
    } else {    # If we are ready to move on to the next Query ID
        if ( $lineage->{tophit} ) {

            # We have finished a hit
            &_print_hits($lineage);
            $finished_id = 1;

            # Initialize new lineage query ID provided our taxon lineage matches what's in the lineage hash
            if ( $tax->{lineage} =~ /$lineage->{lineage}/ ) {
                $lineage->{id}        = $fields->[0];
                $lineage->{best_e}    = $fields->[10];
                $lineage->{best_rows} = [$fields];
            } else {
                $lineage->{best_rows} = [];
                $lineage->{id}        = '';
                $lineage->{best_e}    = 100;
            }
        }
    }

    if ( &_filter_hit( $tax->{name} ) ) {
        push( @$filter_hits, $fields );
    }
    return $finished_id;
}

# Print out the best hit for this particular query ID.
sub _print_hits {
    my $lineage = shift;
    if ( $lineage->{id} ) {
        # We have finished a hit
        if ( scalar @{ $lineage->{best_rows} } ) {
            if ( !&_filter_best_hits( $lineage->{best_rows} ) ) {
                map {
                    print { $lineage->{handle} } join( "\t", @$_ );
                    print { $lineage->{handle} } "\n";
                } @{ $lineage->{best_rows} };
            } else {
                print STDERR "Filtered $lineage->{name}\n";
            }
        }
    }
}

# Append trace information to the m8 hit.
# Adds 'template_id', which is assuned to the query name if trace mapping file doesn't exist
# Adds 'trace_end', which is either forward or reverse.
sub _add_trace_info {
    my $list = shift;

    my $id = $list->[0];
    if ( $trace_lookup->{ $id } ) {
        push( @$list,
            (
                $trace_lookup->{$id}->{'template_id'},
                $trace_lookup->{$id}->{'trace_end'}
            )
        );
    } elsif ( !$trace_mapping_file ) {

        # Ghetto way of checking for directionality.
        my $dir = 'F';
        if ( $id =~ /(.*)[\_\/](\d)/ ) {
			# Determine forward and reverse by which mate pair it is
            if ( $2 == 1 ) {
                push( @$list, ( $1, 'F' ) );
            }
            elsif ( $2 == 2 ) {
                push( @$list, ( $1, 'R' ) );
            }
            else {
                print STDERR "Couldn't figure out the clone name from $list->[0] assuming F\n";
                push( @$list, ( $1, 'F' ) );
            }
        }
    } else {
        print STDERR "couldn't find trace info for $list->[0]\n";
    }
    return $list;
}

# Ensure query and subject have enough overlap to be a great hit
sub _filter_best_hits {
    my $hits   = shift;
    my $filter = 0;
    foreach my $fhit (@$filter_hits) {
        foreach my $hit (@$hits) {
            my $overlap = &_get_overlap_length( [ $fhit->[6], $fhit->[7] ],
                [ $hit->[6], $hit->[7] ] );
            if ( $overlap >= $min_overlap ) {
                print STDERR
                  "Here to filter out $hit->[0] $fhit->[14] with $overlap\n";
                $filter = 1;
                last;
            }
        }
        last if $filter;
    }
    return $filter;

}

# If argument was provided to filter a specific lineage, then do it
sub _filter_hit {
    my $lineage = shift;

    my $retval = 0;
    if ( $filter_lineage && $lineage =~ /$filter_lineage/i ) {
        print STDERR "Filtering out because of $lineage\n";
        $retval = 1;
    }
    return $retval;
}

sub _get_overlap_length {
    my ( $q, $s ) = @_;

    my $qlen = $q->[1] - $q->[0];
    my $slen = $s->[1] - $s->[0];
    my $min  = min( $q->[1], $q->[0], $s->[1], $s->[0] );
    my $max  = max( $q->[1], $q->[0], $s->[1], $s->[0] );

    my $len = ( $qlen + $slen ) - ( $max - $min );

    return $len;
}

sub _create_list_of_outputs {
    my $config = shift;
    my $out_dir;
    my ( $name, $directories, $suffix );
    if ( $config->{output_dir} ) {
		$out_dir = $config->{output_dir};
	}

	if ( $config->{fasta} ) {
        ( $name, $directories, $suffix ) =
          fileparse( $config->{fasta}, qr/\.[^.]*/ );
        $out_dir = $directories if !defined $out_dir;
    } elsif ( $config->{blast} ) {
        ( $name, $directories, $suffix ) =
          fileparse( $config->{blast}, qr/\.[^.]*/ );
        $out_dir = $directories if !defined $out_dir;
    } else {
        confess "Must pass &BestBlast2 an output_dir. $!\n";
    }
    open OUT, ">$out_dir/$name\_filtered_blast.list"
      or confess "Unable to open $out_dir/$name\_filtered_blast.list\n";
    print OUT "$overallfile\n";
    print OUT "$out1file\n";
    print OUT "$out2file\n";
    close OUT;
    return "$out_dir/$name\_filtered_blast.list";
}

1;
