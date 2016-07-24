
=head1 NAME

LGTsam2lca.pm - Run bwa and generate lca's

=head1 SYNOPSIS

Need to put something useful here

=head1 DESCRIPTION

A module to run bwa

=head1 AUTHOR - David R. Riley

e-mail: driley@som.umaryland.edu

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _
my $results = GetOptions (
    \%options,
    'file_list=s',
    'out_file=s',
    'ncbitax=s',
    'gitax=s',
    'dbhost=s',
    'taxondb=s',
    'taxoncoll=s',
    'overwrite=s',
    'idx_dir=s',
    'chunk_size:s',
    'taxon_dir:s',
    'samtools_bin:s',
    'step:s',
    'check_mates:s',
    'independent_lca',
    'se_lca',
    'help|h') || pod2usage();

=cut

package LGT::LGTsam2lca;
use strict;
use warnings;
use version;
use File::Basename;
use LGT::Common;

$| = 1;

my $check_mates;
my $CHUNK_SIZE = 10000;

my $chunk   = [];
my $unmapped_counts = 0;
my $step    = 0;
my $out2;
my $samtools;
my $outf;

my %options;

=head2 new

 Title   : new
 Usage   : my $lgtsam2lca = LGTsam2lca->new({out_file => $outfile,...})
 Function: Creates a new LGTsam2lca object.
 Returns : An instance of LGTsam2lca
 Args    : A hash containing potentially several config options:

           out_file - Path to the output file.
           samtools_bin - path to the samtools binary
           gi2tax - GiTaxon object
=cut

sub new {
    my ( $class, $args ) = @_;

    my $step = $args->{step} ? $args->{step} : 0;
    my $self = {
        reads_by_read_id => {},
        reads_by_mate_id => {},
        out_file         => $args->{out_file},
        samtools_bin     => $args->{samtools_bin},
        gi2tax           => $args->{gi2tax},
        step             => $step,
        check_mates      => $args->{check_mates},
        complete_bam     => $args->{complete_bam}
    };
    bless $self, $class;

    $samtools = $self->{samtools_bin};
    if ( $self->{complete_bam} ) {
        $self->prime_hash( $self->{complete_bam} );
    }
    return $self;
}

=head2 runSam2Lca

 Title   : runSam2Lca
 Usage   : my $output_files = $lgtsam2lca->runSam2lca()
 Function: run sam2lca on the input sam/bam list
 Returns : An object with the output files in it
 Args    :
=cut

sub runSam2Lca {
    my $self   = shift;
    my $config = shift;

    %options = %$config;

    # If the user passes in a list of files, then we'll process them here.
    if ( $self->{file_list} ) {

        my $count = `wc -l $self->{file_list}`;
        $step = $self->{step};
        print STDERR "Working on "
          . ( $step * $CHUNK_SIZE ) . " to "
          . ( $step * $CHUNK_SIZE + $CHUNK_SIZE ) . "\n";
        open IN, "<$self->{file_list}"
          or die "Couldn't open $self->{file_list}\n";
        my $cnter = 0;
        while (<IN>) {
            chomp;
            $self->process_file( { file => $_ } );
        }
        print STDERR
          "\nFinished looping through step $step. Now writing output\n";
        return $self->writeOutput();

    } else {
        print STDERR "Called runSaml2Lca with no file list, nothing to do!\n";
    }
}

=head2 writeOutput

 Title   : writeOutput
 Usage   : my $output_files = $lgtsam2lca->writeOutput()
 Function: Write the output of the currently compiled LCAs
 Returns : An object with the output files in it.
 Args    :
=cut

sub writeOutput {
    my $self = shift;
    my $out2;
    my $outf = $self->{out_file} ? ">$self->{out_file}" : '>-';

    #    if($self->{independent_lca}){
    if ( $self->{out_file} ) {
        $outf =~ /^>(.+?)\.(\w+?)$/;
        $out2 = "$1\_independent_lca.$2";
    } else {
        $out2 = ">-";
    }
    open( OUT2, ">", "$out2" ) or die "Couldn't open $out2\n";

    #    }

    open my $out, $outf or die "Couldn't open output\n";
    foreach my $key ( keys %{ $self->{reads_by_mate_id} } ) {
        #print STDERR "Mate seen:$key processing with ...";

		# Printing out the independent LCA for each read
        #   print STDERR "\tindependent_lca ...";}
		if ( !defined( $self->{reads_by_read_id}->{"$key\_1"} ) ) {
			$self->{reads_by_read_id}->{"$key\_1"} = '';
		}
        print OUT2
          join( "\t", ( "$key\_1", $self->{reads_by_read_id}->{"$key\_1"} ) );
        print OUT2 "\n";

        if ( !defined( $self->{reads_by_read_id}->{"$key\_2"} ) ) {
			$self->{reads_by_read_id}->{"$key\_2"} = '';
		}
        print OUT2
          join( "\t", ( "$key\_2", $self->{reads_by_read_id}->{"$key\_2"} ) );
        print OUT2 "\n";

		# If mate has no LCA it didn't map to the reference genome
        if ( !defined( $self->{reads_by_mate_id}->{$key} ) ) {
			#print STDERR "Found no LCA for mate $key ... skipping\n";
			$self->{reads_by_mate_id}->{$key} = '';
			#next;
        }

		# Determine the conservative single-end LCA and write to conservative outfile
        #      print STDERR "\tSingleEnd_lca ...";
		my $new_conservative_se_lca = &find_lca(
            [
                $self->{reads_by_read_id}->{"$key\_1"},
                $self->{reads_by_read_id}->{"$key\_2"}
            ]
        );

        if ( $self->{reads_by_read_id}->{"$key\_1"} =~
               $self->{reads_by_read_id}->{"$key\_2"}
            || $self->{reads_by_read_id}->{"$key\_2"} =~
            $self->{reads_by_read_id}->{"$key\_1"} )
        {
            ## Liberal
            if (
                length( $self->{reads_by_read_id}->{"$key\_1"} ) >=
                length( $self->{reads_by_read_id}->{"$key\_2"} ) )
            {
                print $out join(
                    "\t",
                    (
                        $key,
                        $self->{reads_by_mate_id}->{$key},
                        $new_conservative_se_lca,
                        $self->{reads_by_read_id}->{"$key\_1"}
                    )
                );
                print $out "\n";
            } else {
                print $out join(
                    "\t",
                    (
                        $key,
                        $self->{reads_by_mate_id}->{$key},
                        $new_conservative_se_lca,
                        $self->{reads_by_read_id}->{"$key\_2"}
                    )
                );
                print $out "\n";
            }
        } else {
            print $out join(
                "\t",
                (
                    $key,
                    $self->{reads_by_mate_id}->{$key},
                    $new_conservative_se_lca,
                    $new_conservative_se_lca
                )
            );
            print $out "\n";
        }
    }

    return {
        independent => $out2,
        normal      => $outf
    };
}

sub prime_hash {
    my $self = shift;
    my $f    = shift;
    print STDERR
      "LGT::LGTsam2lca --- Complete BAM file provided.  Now priming 'reads by mate id' hash...\n";
    open( my $in, "-|", "$samtools view $f" )
      or die "Unable to open $f for priming the hash\n";
    while (<$in>) {
        my @fields = split(/\t/);
        $self->{reads_by_mate_id}->{ $fields[0] } = undef;
    }
    print STDERR "LGT::LGTsam2lca --- Finished priming hash\n";
}

sub process_sam_line {
    my ( $self, $l ) = @_;
    chomp $l;

    # Don't count @seq lines
    if ( $l =~ /^@/ ) {
        next;
    }

    #$count++;
    my @fields    = split( /\t/, $l );
    my $flag      = parse_flag( $fields[1] );
    my $read_name = "$fields[0]\_1";
    if ( !$flag->{first} ) {
        $read_name = "$fields[0]\_2";
    }

    # Keep track of how many mate pairs had both mates unmapped
	$unmapped_counts++ if ($flag->{qunmapped} && $flag->{munmapped} && $flag->{first});

	# If both mates fail to map to the reference don't add to hash
	return if ($flag->{qunmapped} && $flag->{munmapped});

	# Here we determine LCA for each read ID (2 per mate pair) and for each mate ID

	# If the current read is mapped add to read_id hash
    if(!$flag->{qunmapped}) {
		my $tax = $self->{gi2tax}->getTaxon( $fields[2] );
		#print STDERR "No lineage found for $fields[2]\n" if (! defined $tax->{lineage});

		# Here we'll deal with keeping track of things by read
	    if ( $self->{reads_by_read_id}->{$read_name} ) {
	        my $lca = &find_lca(
	            [ $self->{reads_by_read_id}->{$read_name}, $tax->{lineage} ] );
	        $self->{reads_by_read_id}->{$read_name} = $lca;
	        #print "$read_name\t$tax->{lineage}\t$lca\n";
	    } else {
	        $self->{reads_by_read_id}->{$read_name} = $tax->{lineage};
	    }

	    # If we a) aren't checking mates or b) are and both mates map to reference...
	    if (  !$self->{check_mates}
	        || ($self->{check_mates} && !$flag->{munmapped}) )
	    {
	        # Here we'll keep track of things by mate
	        if ( $self->{reads_by_mate_id}->{ $fields[0] } ) {
	            #print STDERR "$fields[0] are found\n";
	            my $lca = &find_lca(
	                [
	                    $self->{reads_by_mate_id}->{ $fields[0] },
	                    $tax->{lineage}
	                ]
	            );
	            $self->{reads_by_mate_id}->{ $fields[0] } = $lca;
	        } else {
	            #print STDERR "$fields[0] are not found\n";
	            $self->{reads_by_mate_id}->{ $fields[0] } = $tax->{lineage};
	        }
		}
	}
}

sub process_file {
    my ( $self, $config ) = @_;
    my $line = $config->{file};
    my $handle;
    my $presplit = 0;
    my $file     = $line;

    if ( $config->{handle} ) {
        $handle = $config->{handle};
        print STDERR "LGT::LGTsam2lca" . " ---File handle provided\n";
        open $handle, "<$file" or die "Unable to open $file\n";
    } elsif ( $file =~ /.bam$/ ) {
        print STDERR "LGT::LGTsam2lca" . " ---BAM file provided\n";
        open( $handle, "-|", "$samtools view $file" )
          or die "Unable to open $file\n";
    } else {
        die "Need to pass a valid file or a valid filehandle\n";
    }

	print STDERR "LGT::LGTsam2lca --- now processing BAM file $file\n";
    # Loop till we're done.
    my $end   = $CHUNK_SIZE;
    my $count = 0;
    my $l;
    my $hit = 0;

    while ( $l = <$handle> ) {
        chomp $l;
        # Don't count @seq lines
        if ( $l =~ /^@/ ) {
            next;
        }
        $count++;
        $self->process_sam_line($l);
    }
	print STDERR "There were $unmapped_counts reads that did not map to the reference genome\n";
	print STDERR "LGT::LGTsam2lca --- finished processing BAM file $file\n";
}

1;
