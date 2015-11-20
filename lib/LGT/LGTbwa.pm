
=head1 NAME

LGTbwa - Run bwa.

=head1 SYNOPSIS

Need to put something useful here

=head1 DESCRIPTION

A module to run bwa

=head1 AUTHOR - David R. Riley

e-mail: driley@som.umaryland.edu

=head1 APPENDIX

The rest of the documentation details each of the object methods.
Internal methods are usually preceded with a _

$results = GetOptions (\%options, 
              'input_dir|i:s',
              'input_base|ib:s',
              'input_bam:s',
              'paired=i',
              'ref_file|r=s',
              'ref_file_list|rl=s',
              'output_dir|o=s',
              'output_bam:s',
              'samtools_flag:s',
              'samtools_path:s',
              'mismatch|mm=s',
              'max_gaps|mg=s',
              'max_gap_extentions|mge=s',
              'open_gap_penalty|og=s',
              'extend_gap_penalty|eg=s',
              'threads|t=s',
              'use_bwasw',
              'num_aligns|na=s',
              'bwa_path|b=s',
              'cleanup_sai:s',
              'run_lca:s',
              'samtools_opts:',
              'help|h');

=cut

package LGT::LGTbwa;
use strict;
use warnings;
use Pod::Usage;
use File::Basename;
use File::Find;
my %options   = ();
my $ref_files = [];
my $PAIRED;
my $options_string = '';
my $sam2lca;

sub runBWA {
    my $opts = shift;
    %options = %$opts;

    ## make sure all passed options are peachy
    &check_parameters( \%options );

    # This is a HACK right now.
    if ( $options{run_lca} ) {
        if ( -e $options{out_file} && $options{overwrite} == 0 ) {
            print STDERR
              "Already found the output for BWA-LCA: $options{out_file}\n";
            return $options{out_file};
        }
        my $lca_conf = {
            samtools_bin => "$options{samtools_path}/samtools",
            out_file     => $options{out_file},
            gi2tax       => $options{lgtseek}->getGiTaxon( {} )
        };
        if ( $options{input_bam} ) {
            $lca_conf->{complete_bam} = $options{input_bam};
        }
        $sam2lca = LGT::LGTsam2lca->new($lca_conf);
    }

    $PAIRED = defined $options{paired} ? $options{paired} : 1;

    if ( $options{use_bwasw} ) {
        die "bwasw is not implemented yet\n";
    }
    else {
        &run_bwa();
    }
}

sub run_bwa {

    foreach my $ref (@$ref_files) {
        chomp $ref;
        my $refname = '';
        if ( $ref =~ /.*\/([^\/]+)\.[^\.]+$/ ) {
            $refname = $1;
        }
        else {
            $ref =~ /.*\/([^\/]+)\.?[^\.]*$/;
            $refname = $1;
        }

        if ( !-e "$ref.bwt" ) {
            print STDERR "Did not find reference file $ref.bwt\n";

            print STDERR "Indexing $ref\n";
            my $cmd = "$options{bwa_path} index $ref";
            system($cmd) == 0 or die "Unable to index $ref\n";
        }

        # In here if we're paired
        if ($PAIRED) {
            my $in1;
            my $in2;
            my $out1;
            my $out2;

            open( my $bwa_log, ">>$options{output_dir}/\.bwa_stderr.log" )
              || die
              "&run_bwa can't open: $options{output_dir}/\.bwa_stderr.log";
            my $output_file =
              "$options{output_dir}/$refname\_$options{input_base}.sam";
            if ( $options{output_bam} ) {
                $output_file =
                  "$options{output_dir}/$refname\_$options{input_base}.bam";
            }

            if ( ( -e $output_file ) && !$options{overwrite} ) {
                print STDERR
                  "Already found the output for &runBWA: $output_file.\n";
                return;
            }
            if ( $options{input_bam} ) {
                my ( $name, $dir, $suff ) =
                  fileparse( $options{input_bam}, ( "_prelim.bam", ".bam" ) );
                $in1                 = $options{input_bam};
                $in2                 = $options{input_bam};
                $options{input_base} = $name;

                # Run the first one through aln
                $out1 = "$options{output_dir}/$refname\_$name\_aln_sa1.sai";
                my $cmd =
"$options{bwa_path} aln -b -1 $options_string $ref $options{input_bam} > $out1 2>>$options{output_dir}/\.bwa_stderr.log";
                print $bwa_log "CMD: $cmd\n";
                system($cmd) == 0 or die "Unable to run $cmd\n";

                # Run the second one through aln
                $out2 = "$options{output_dir}/$refname\_$name\_aln_sa2.sai";
                $cmd =
"$options{bwa_path} aln -b -2 $options_string $ref $options{input_bam} > $out2 2>>$options{output_dir}/\.bwa_stderr.log";
                print $bwa_log "CMD: $cmd\n";
                system($cmd) == 0 or die "Unable to run $cmd\n";
            }
            else {
                find(
                    sub {
                        if ( $File::Find::name =~
/$options{input_base}(_[12]{1}){0,1}(\.\w+)?\.f\w{0,3}q(.gz)?/
                          )
                        {
                            $in1 = $File::Find::name if ( $1 =~ /_1/ );
                            $in2 = $File::Find::name if ( $1 =~ /_2/ );
                        }
                    },
                    ( $options{input_dir} )
                );

                $out1 =
"$options{output_dir}/$refname\_$options{input_base}_aln_sa1.sai";
                $out2 =
"$options{output_dir}/$refname\_$options{input_base}_aln_sa2.sai";

                # Run the first one through aln
                if ( $options{overwrite} || !-e $out1 ) {
                    my $cmd =
"$options{bwa_path} aln $options_string $ref $in1 > $out1 2>>$options{output_dir}/\.bwa_stderr.log";
                    print $bwa_log "CMD: $cmd\n";
                    system($cmd) == 0 or die "Unable to run $cmd\n";
                }
                if ( $options{overwrite} || !-e $out2 ) {

                    # Run the second one through aln
                    my $cmd =
"$options{bwa_path} aln $options_string $ref $in2 > $out2 2>>$options{output_dir}/\.bwa_stderr.log";
                    print $bwa_log "CMD: $cmd\n";
                    system($cmd) == 0 or die "Unable to run $cmd\n";
                }
            }

            my $cmd;

            # Run sampe
            if ( $options{run_lca} ) {
                my $cmd =
"$options{bwa_path} sampe -n $options{num_aligns} $ref \"$out1\" \"$out2\" \"$in1\" \"$in2\" 2>>$options{output_dir}/\.bwa_stderr.log | $options{samtools_path}samtools view $options{samtools_flag} -S - 2>>$options{output_dir}/\.bwa_stderr.log";
                open( my $handle, "-|", $cmd );
                $sam2lca->process_file( { handle => $handle } );
            }
            elsif ( $options{output_bam} ) {
                if ( $options{overwrite}
                    || !
                    -e "$options{output_dir}/$refname\_$options{input_base}.bam"
                  )
                {
                    my $cmd =
"$options{bwa_path} sampe -n $options{num_aligns} $ref \"$out1\" \"$out2\" \"$in1\" \"$in2\" 2>>$options{output_dir}/\.bwa_stderr.log | $options{samtools_path}samtools view $options{samtools_flag} -bS - > $options{output_dir}/$refname\_$options{input_base}.bam 2>>$options{output_dir}/\.bwa_stderr.log";
                    print $bwa_log "CMD: $cmd\n";
                    system($cmd) == 0 or die "Unable to run $cmd\n";
                }
            }
            else {
                if ( $options{overwrite}
                    || !
                    -e "$options{output_dir}/$refname\_$options{input_base}.sam"
                  )
                {
                    my $cmd =
"$options{bwa_path} sampe -n $options{num_aligns} $ref \"$out1\" \"$out2\" \"$in1\" \"$in2\" > $options{output_dir}/$refname\_$options{input_base}.sam 2>>$options{output_dir}/\.bwa_stderr.log";
                    print $bwa_log "CMD: $cmd\n";
                    system($cmd) == 0 or die "Unable to run $cmd\n";
                }
            }

            if ( $options{cleanup_sai} ) {
                my $cmd = "rm -f $out1 $out2";
                print $bwa_log "CMD: $cmd\n";
                system($cmd) == 0 or die "Unable to run $cmd\n";
            }
        }

        # In here if we aren't paired
        else {
            my $single_read_file =
              "$options{output_dir}/$options{input_base}_single_read.txt";
            my $cmd = "touch $single_read_file";
            system($cmd) == 0 or die "Unable to run $cmd\n";

            my $in;
            find(
                sub {
                    $in = $File::Find::name
                      if ( $File::Find::name =~
/$options{input_base}(_[12]{1})?\.\w+)?\.f\w{0,3}q(.gz)?/
                      );
                },
                ( $options{input_dir} )
            );
            my $out =
              "$options{output_dir}/$refname\_$options{input_base}_aln_sa.sai";
            $cmd = "$options{bwa_path} aln $options_string $ref $in > $out";
            if ( $options{overwrite} || !-e "$out" ) {
                system($cmd) == 0 or die "Unable to run $cmd\n";
            }
            if ( $options{overwrite}
                || !
                -e "$options{output_dir}/$refname\_$options{input_base}.sam" )
            {
                $cmd =
"$options{bwa_path} samse -n $options{num_aligns} $ref \"$out\" \"$in\" > $options{output_dir}/$refname\_$options{input_base}.sam";
                system($cmd) == 0 or die "Unable to run $cmd\n";
            }
            if ( $options{cleanup_sai} ) {
                $cmd = "rm -f $out";
                system($cmd) == 0 or die "Unable to run $cmd\n";
            }
        }
    }

    if ( $options{run_lca} ) {
        $sam2lca->writeOutput();
    }
}

sub check_parameters {

    ## make sure input file exists
    if ( !-e $options{'input_dir'} && !-e $options{'input_bam'} ) {
        print STDERR "Input invalid\n";
        pod2usage( { -exitval => 0, -verbose => 0, -output => \*STDOUT } );
    }

    if ( $options{ref_file} ) {
        my @files = split( /,/, $options{ref_file} );
        $ref_files = \@files;
    }
    elsif ( $options{ref_file_list} ) {
        my @lines = `cat $options{ref_file_list}`;
        $ref_files = \@lines;
    }
    else {
        die "No reference file specified in ref_file or ref_file_list";
    }

    my $options_to_param_keys = {
        'mismatch'           => '-M',
        'max_gaps'           => '-o',
        'max_gap_extensions' => '-e',
        'open_gap_penalty'   => '-O',
        'extend_gap_penalty' => '-E',
        'threads'            => '-t'
    };

    my $opts = [];
    foreach my $key ( keys %$options_to_param_keys ) {
        if ( $options{$key} ) {
            push( @$opts, "$options_to_param_keys->{$key} $options{$key}" );
        }
    }
    $options_string = join( " ", @$opts );

    return 1;
}
1;
