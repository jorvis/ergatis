
=head1 NAME

LGTFinder - Find putative LGT in BLAST data

=head1 SYNOPSIS

Need to put something useful here

=head1 DESCRIPTION

A module to take output from LGTBestBlast and look for putative LGT's both
within reads and across read pairs

=head1 AUTHOR - David R. Riley

e-mail: driley@som.umaryland.edu

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

=cut

package LGT::LGTFinder;

use strict;
use warnings;
use LGT::Common;
use File::Basename;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use List::Util qw[min max];
$| = 1;

my $traces_by_template = {};
my $traces_by_trace_id = {};

my $MAX_OVERLAP     = 20;
my $MIN_LENGTH      = 0;
my $DEFAULT_LINEAGE = 'Homo';	# This is only applicable for human/bacteria LGT

my $output_dir;
my $filename;

my ( $max_overlap_threshold, $min_len_threshold, $ref_lineage );

# Initial script to process input, setup output filehandles, and call the find_lgt function
sub findLGT {
    my ( $class, $options ) = @_;

    $max_overlap_threshold =
      defined $options->{max_overlap} ? $options->{max_overlap} : $MAX_OVERLAP;
    $min_len_threshold =
      defined $options->{min_length} ? $options->{min_length} : $MIN_LENGTH;
    $ref_lineage =
      defined $options->{ref_lineage}
      ? $options->{ref_lineage}
      : $DEFAULT_LINEAGE;

    ## Determine output_dir and prefix
    $filename = 'lgt_finder';
    if ( $options->{input_file_list} ) {
        my ( $name, $directories, $suffix ) =
          fileparse( $options->{input_file_list}, qr/\.[^.]*/ );
        $output_dir =
          defined $options->{output_dir}
          ? $options->{output_dir}
          : $directories;
        $filename = $name;
    } elsif ( $options->{input} ) {
        my ( $name, $directories, $suffix ) =
          fileparse( $options->{input}, qr/\.[^.]*/ );
        $output_dir =
          defined $options->{output_dir}
          ? $options->{output_dir}
          : $directories;
    }
    if ( $options->{output_prefix} ) {
        $filename = $options->{output_prefix};
    }

    # Open up our filehandles and write the headers
    open OUTCLONE, ">$output_dir/$filename\_by_clone.txt"
      or die
      "Couldn't open by_clone output: $output_dir/$filename\_by_clone.txt $!";
    open OUTTRACE, ">$output_dir/$filename\_by_trace.txt"
      or die
      "Couldn't open by_trace output: $output_dir/$filename\_by_trace.txt $!";
    print OUTTRACE
      "Read\t\tDonor-Eval\tLength\tStart\tend\tRef-Eval\tLength\tStart\tend\tDonor-Lineage\tRef-Lineage\n";

    # Handle differently single input file vs input file list
    if ( $options->{input} ) {
        my $old_prefix;
        foreach my $file ( split( /,/, $options->{input} ) ) {
            chomp;

            # Get prefix name
            my ( $fn, $directories, $suffix ) = fileparse( $file, qr/\.[^.]*/ );
            my $prefix = $options->{output_prefix};
            if ( !$prefix ) {
                $fn =~ /^(.*)\_[^\_]+/;
                $prefix = $1;
            }

           # Find LGT in donor or host file just before processing the next file
            if ( $prefix && ( $prefix ne $old_prefix ) ) {
                &_find_lgt_in_trace();
                $traces_by_trace_id = {};
                $old_prefix         = $prefix;
            }
            print STDERR "======== &LGTFINDER: Processing $file ========\n";
            &_process_file($file, $options);
        }

        # Get LGT in clones, and in the final trace
        &find_lgt();
        close OUTCLONE;
        close OUTTRACE;
    } elsif ( $options->{input_file_list} ) {
        open IN, "<$options->{input_file_list}"
          or die "Unable to open input $options->{input_file_list}\n";
        my ( $filename, $directories, $suffix ) =
          fileparse( $options->{input_file_list}, qr/\.[^.]*/ );

        my $old_prefix;
        while (<IN>) {
            chomp;
            my ( $filename, $directories, $suffix ) =
              fileparse( $_, qr/\.[^.]*/ );
            $filename =~ /^([^.]*)_[^_]+$/;
            my $prefix = $1;

           # Find LGT in donor or host file just before processing the next file
            if ( $prefix && ( $prefix ne $old_prefix ) ) {
                &_find_lgt_in_trace();
                $traces_by_trace_id = {};
                $old_prefix         = $prefix;
            }
            print STDERR "======== &LGTFINDER: Processing $_ ========\n";
            &_process_file($_, $options);
        }

        # Get LGT in clones (overall) and in the final trace (donor/host)
        &find_lgt();
        close OUTCLONE;
        close OUTTRACE;
    }

    my $valid_count =
      `grep ';$options->{lineage1};' $output_dir/$filename\_by_clone.txt | grep ';$options->{lineage2};' | wc -l`;
    my $valid_int_count =
      `wc -l $output_dir/$filename\_by_trace.txt | cut -f1 -d \' \'`;
    chomp $valid_count;
    chomp $valid_int_count;
    return {
        valid_clones => $valid_count,
        valid_traces => $valid_int_count,
        by_clone     => "$output_dir/$filename\_by_clone.txt",
        by_trace     => "$output_dir/$filename\_by_trace.txt",
    };
}

sub _process_file {
    my $file = shift;
	my $opts = shift;
    $file =~ /_([^_]+).out$/;
    my $type = $1;
    if ( $type eq 'overall' ) {
        &_process_overall($file, $opts);
    } else {
        &_process_within_trace( $file, $type );
    }
}

# Gather trace information for each template ID
sub _process_overall {
    my $file = shift;
	my $options = shift;
    open IN2, "<$file" or die "Unable to open $file\n";
	# Read in hits, line by line
    while (<IN2>) {
        chomp;
        my @fields = split( /\t/, $_ );

        my $hit = {
            trace_id        => $fields[0],
            accession       => $fields[1],
            pid             => $fields[2],
            align_len       => $fields[3],
            mismatches      => $fields[4],
            gaps            => $fields[5],
            query_start     => $fields[6],
            query_end       => $fields[7],
            subj_start      => $fields[8],
            subj_end        => $fields[9],
            evalue          => $fields[10],
            score           => $fields[11],
            taxon_id        => $fields[12],
            scientific_name => $fields[13],
            lineage         => $fields[14],
            template_id     => $fields[15],		# may be query_id (trace_id)
            direction       => $fields[16]
        };

        my $dir = 'forward';
        if ( $hit->{direction} eq 'R' ) {
            $dir = 'reverse';
        }

        $hit->{scientific_name} =~ /^(\w+) /;
        my $genera = $1;

        # Start initializing trace hash for given template ID if it doesn't exist
        if ( !$traces_by_template->{ $hit->{template_id} } ) {
            $traces_by_template->{ $hit->{template_id} } = {
                'forward' => {},
                'reverse' => {}
            };
        }
        if ( !$traces_by_template->{ $hit->{template_id} }->{$dir}->{$genera} ) {
            $traces_by_template->{ $hit->{template_id} }->{$dir}->{$genera} =[];
        }
		
		# Filter out certain genus names and process the remaining
        if ( !&filter_genera($genera) ) {
            $traces_by_template->{ $hit->{template_id} }->{$dir}->{traces}
              ->{ $hit->{trace_id} } = 1;
            $hit->{genera} = $genera;

            my $hit_filter =
              $traces_by_template->{ $hit->{template_id} }->{$dir}->{hit}
              ->{hit_filter}
              ? $traces_by_template->{ $hit->{template_id} }->{$dir}->{hit}
              ->{hit_filter}
              : 0;

            my $lca =
              defined(
                $traces_by_template->{ $hit->{template_id} }->{$dir}->{hit}->{lca} )
              ? $traces_by_template->{ $hit->{template_id} }->{$dir}->{hit}->{lca}
              : $hit->{lineage};

			# Mark for hit filtering if lineage/tag options were provided
            if (   $options->{lca_filter_lineage}
                && $options->{lca_filter_tag}
                && $hit->{lineage} =~ /$options->{lca_filter_lineage}/
                && $hit->{trace_id} =~ /$options->{lca_filter_tag}/ )
            {
                $hit_filter = 1;
                $lca =
                  defined(
                    $traces_by_template->{ $hit->{template_id} }->{$dir}->{hit}->{lca} )
                  ? $traces_by_template->{ $hit->{template_id} }->{$dir}->{hit}->{lca}
                  : undef;
            } else {
                if ( !$traces_by_template->{ $hit->{template_id} }->{$dir}
                    ->{lineages}->{ $hit->{lineage} } )
                {
                    $traces_by_template->{ $hit->{template_id} }->{$dir}
                      ->{lineages}->{ $hit->{lineage} } = 0;
                }
                $traces_by_template->{ $hit->{template_id} }->{$dir}
                  ->{lineages}->{ $hit->{lineage} }++;
                $lca = &find_lca( [ $hit->{lineage}, $lca ] );
            }

			# Once hit_filter and LCA are determined, fill in missing info
            $traces_by_template->{ $hit->{template_id} }->{$dir}->{hit} = {
                evalue     => $hit->{evalue},
                align_len  => $hit->{align_len},
                trace_id   => $hit->{trace_id},
                genera     => $hit->{genera},
                hit_filter => $hit_filter,
                lca        => $lca
            };

            push(
                @{
                    $traces_by_template->{ $hit->{template_id} }->{$dir}
                      ->{genera}->{$genera}
                },
                $hit->{trace_id}
            );
        }
    }
    close IN2;
}

sub find_lgt {
    if ( keys %$traces_by_template ) {
        print STDERR
          "======== &LGTFINDER: Going to look for LGT within clones/mates\n";
        &_find_lgt_in_clone();
    }
    if ( keys %$traces_by_trace_id ) {
        &_find_lgt_in_trace();
    }
}

sub filter_genera {
    my $genera = shift;

    my $genera_to_ignore = {
        'Pan'        => 1,
        'Papio'      => 1,
        'Macaca'     => 1,
        'synthetic'  => 1,
        'Synthetic'  => 1,
        'uncultured' => 1,
        'Uncultured' => 1
    };

    my $retval = 0;
    if ( $genera_to_ignore->{$genera} ) {
        $retval = 1;
    }
    return $retval;
}

sub _find_lgt_in_trace {

    foreach my $trace ( keys %$traces_by_trace_id ) {
        my $found_lgt     = 0;
        my $genera        = {};
        my $evalue        = '';
        my $length        = '';
        my $ref_hit;
        my $donor_hit;

        # Look at all combinations of genera to see if any have LGT evidence.
        foreach my $g1 ( keys %{ $traces_by_trace_id->{$trace} } ) {
            foreach my $h1 ( @{ $traces_by_trace_id->{$trace}->{$g1} } ) {
                foreach my $g2 ( keys %{ $traces_by_trace_id->{$trace} } ) {
                    foreach my $h2 ( @{ $traces_by_trace_id->{$trace}->{$g2} } )
                    {
                        if ( ( $g1 ne $g2 )
                            && &check_overlap( $h1, $h2 )
                            && ( $h1->{align_len} >= $min_len_threshold )
                            && ( $h2->{align_len} >= $min_len_threshold ) )
                        {
							# Determin which hit is the donor and which is the host
                            if ( $h1->{lineage} =~ /$ref_lineage/ ) {
                                $genera->{"$g2:$g1"} = 1;
                                $found_lgt           = "$g2\t$g1";
                                $evalue              = $h2->{evalue};
                                $length              = $h2->{align_len};
                                $ref_hit             = $h1;
                                $donor_hit           = $h2;
                            } else {
                                $genera->{"$g1:$g2"} = 1;
                                $found_lgt           = "$g1\t$g2";
                                $evalue              = $h1->{evalue};
                                $length              = $h1->{align_len};
                                $ref_hit             = $h2;
                                $donor_hit           = $h1;
                            }
                        }
                    }
                }
            }
        }
		
		# If LGT was found then print all sorts of nice info about it.
        if ($found_lgt) {
#            print STDERR "Found LGT in trace $trace between ".join(',',keys %$genera)."\n";
            print OUTTRACE join(
                "\t",
                (
                    $trace,                  $donor_hit->{evalue},
                    $donor_hit->{align_len}, $donor_hit->{query_start},
                    $donor_hit->{query_end}, $ref_hit->{evalue},
                    $ref_hit->{align_len},   $ref_hit->{query_start},
                    $ref_hit->{query_end},   $donor_hit->{lineage},
                    $ref_hit->{lineage}
                )
            );
            print OUTTRACE "\n";
        }
        $found_lgt = 0;
    }
}

sub check_overlap {
    my ( $hit1, $hit2 ) = @_;
    my $overlap = &get_overlap_length(
        [ $hit1->{query_start}, $hit1->{query_end} ],
        [ $hit2->{query_start}, $hit2->{query_end} ]
    );

#    print STDERR "here with $hit1->{query_start} $hit1->{query_end} $hit2->{query_start} $hit2->{query_end} $overlap\n";

    my $retval = 0;
    if ( $overlap < $max_overlap_threshold ) {
        $retval = 1;
    }
    return $retval;

}

sub get_overlap_length {
    my ( $q, $s ) = @_;

    my $qlen = $q->[1] - $q->[0];
    my $slen = $s->[1] - $s->[0];
    my $min  = min( $q->[1], $q->[0], $s->[1], $s->[0] );
    my $max  = max( $q->[1], $q->[0], $s->[1], $s->[0] );

    my $len = ( $qlen + $slen ) - ( $max - $min );

    return $len;
}

sub _find_lgt_in_clone {
    my $lgt_by_genera = {};
    foreach my $clone ( keys %$traces_by_template ) {
        my $forward = join( ",",
            keys %{ $traces_by_template->{$clone}->{forward}->{genera} } );
        my $flineage = join( ",",
            keys %{ $traces_by_template->{$clone}->{forward}->{lineages} } );
        my $reverse = join( ",",
            keys %{ $traces_by_template->{$clone}->{reverse}->{genera} } );
        my $rlineage = join( ",",
            keys %{ $traces_by_template->{$clone}->{reverse}->{lineages} } );

        if ( $forward && $reverse ) {
			
        	my $fwd_traces = [];
        	my $rev_traces = [];

			# Get forward hit information
            my $fhit = $traces_by_template->{$clone}->{forward}->{hit};
            my $nft = scalar keys %{ $traces_by_template->{$clone}->{forward}->{traces} };
            push(
                @$fwd_traces,
                join(
                    "\t",
                    (
                        $fhit->{trace_id},  $fhit->{evalue},
                        $fhit->{align_len}, $fhit->{lca},
                        $fhit->{hit_filter}
                    )
                )
            );

			#template_id	F	trace_id	genera	trace_id,eval,align_len,lca,hit_filter
			#my @Ffields = ( $clone, 'F', $nft, $forward, join( ",", @$fwd_traces ) );
			#my $fline = join( "\t", @Ffields );
			
			# Now for the reverse hit information
            my $rhit = $traces_by_template->{$clone}->{reverse}->{hit};
            my $nrt = scalar keys %{ $traces_by_template->{$clone}->{forward}->{traces} };
            push(
                @$rev_traces,
                join(
                    "\t",
                    (
                        $rhit->{trace_id},  $rhit->{evalue},
                        $rhit->{align_len}, $rhit->{lca},
                        $rhit->{hit_filter}
                    )
                )
            );
			#my @Rfields = ( $clone, 'R', $nft, $reverse, join( ",", @$rev_traces ) );
			#my $rline = join( "\t", @Rfields );

            my $output_line = 0;
			# If our forward hit belongs to the host
            if (   $flineage =~ /$ref_lineage/
                && $rlineage !~ /$ref_lineage/ )
            {
                my @outfields = (
                    $clone, $reverse, $forward,
                    join( ",", @$rev_traces ),
                    join( ",", @$fwd_traces )
                );
                $output_line = join( "\t", @outfields );

                # print STDERR "$fline\n$rline\n";
                push( @{ $lgt_by_genera->{$reverse} }, $clone );
            }
			# If our reverse hit belongs to the host
            elsif ($rlineage =~ /$ref_lineage/
                && $flineage !~ /$ref_lineage/ )
            {
                # print STDERR"$fline\n$rline\n";
                my @outfields = (
                    $clone, $forward, $reverse,
                    join( ",", @$fwd_traces ),
                    join( ",", @$rev_traces )
                );
                $output_line = join( "\t", @outfields );
                push( @{ $lgt_by_genera->{$forward} }, $clone );
            }
			# If neither hit belongs to the host
            else {
                # print "$fline\n$rline\n";
                my @outfields = (
                    $clone, $forward, $reverse,
                    join( ",", @$fwd_traces ),
                    join( ",", @$rev_traces )
                );
                $output_line = join( "\t", @outfields );
                push(
                    @{ $lgt_by_genera->{ $forward . ":" . $reverse } },
                    $clone
                );

            }
            if ($output_line) {
                print OUTCLONE $output_line;
                print OUTCLONE "\n";
            }
        } elsif ($clone) {
            print STDERR "Didn't have both forward and reverse for $clone\n";
        }
    }
}

# Process hits to the donor or host
sub _process_within_trace {
    my $file = shift;
    my $type = shift;
    open IN3, "<$file" or die "Unable to open $file\n";
	# Read hits line by line
    while (<IN3>) {
        chomp;
        my @fields = split( /\t/, $_ );

        my $hit = {
            trace_id        => $fields[0],
            accession       => $fields[1],
            pid             => $fields[2],
            align_len       => $fields[3],
            mismatches      => $fields[4],
            gaps            => $fields[5],
            query_start     => $fields[6],
            query_end       => $fields[7],
            subj_start      => $fields[8],
            subj_end        => $fields[9],
            evalue          => $fields[10],
            score           => $fields[11],
            taxon_id        => $fields[12],
            scientific_name => $fields[13],
            lineage         => $fields[14],
            template_id     => $fields[15],
            direction       => $fields[16]
        };
        $hit->{scientific_name} =~ /^(\w+) /;
        my $genera = $1;

        if ( !$traces_by_trace_id->{ $hit->{trace_id} } ) {
            $traces_by_trace_id->{ $hit->{trace_id} } = {};
        }
        push( @{ $traces_by_trace_id->{ $hit->{trace_id} }->{$genera} }, $hit );
    }
    close IN3;

}

1;
