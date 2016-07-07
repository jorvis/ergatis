#!/usr/local/bin/perl -w

#########################################################################################
#											#
# Name	      : lgt_bwa.pl								#
# Version     : 1.0									#
# Project     : LGT Seek Pipeline							#
# Description : Script to align reads to a reference using BWA				#
# Author      : Sonia Agrawal								#
# Date        : September 1, 2015							#
#											#
#########################################################################################

use strict;

# Extended processing of command line options
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);

# Prints usage message from embedded pod documentation
use Pod::Usage;

# Include packages here
use File::Basename;

#############
# CONSTANTS #
#############

###########
# GLOBALS #
###########
my %options = ();
my @ref_files;

# Log file handle;
my $log_fh;
my ( $ERROR, $WARN, $DEBUG ) = ( 1, 2, 3 );
my %type = ();
my ( $file_count, $exit_code );
my ( $in, $out, $out1, $out2, $file_base, $file_dir, $file_ext, $cmd );
my $fastq_paired = 0;    # Determines if fastq seqs are paired-end or not

################
# MAIN PROGRAM #
################
GetOptions(
    \%options,             'reference|r=s',
    'reference_list|R=s',  'use_mem|m=i',
    'bam_paired|p=i',      'input_dir|I=s',
    'input_file|i=s',      'output_dir|d=s',
    'bwa_path|b=s',        'samtools_path|s=s',
    'misMsc=i',            'maxGapO=i',
    'maxGapE=i',           'gapOsc=i',
    'gapEsc=i',            'nThrds|t=i',
    'maxOcc=i',            'bwa_params|B=s',
    'samtools_params|S=s', 'keep_sam|k=i',
    'log|l=s',             'help|h'
) or pod2usage();

pod2usage( { -exitval => 0, -verbose => 2, -output => \*STDERR } )
  if ( $options{'help'} );

check_parameters( \%options );

#Go though the directory and determine the type of files that are to be aligned
determine_format( \%options, \%type );

check_for_single_sample( \%type );

# Create a list file to store all BAM output
my $bam_list = $options{'output_dir'} . "/bam.list";
open OUT_LIST, ">$bam_list" or print_log_msg($ERROR, "Cannot open $bam_list for writing: $!");

$file_count = keys %type;

foreach my $ref (@ref_files) {
	chomp $ref;
    my $refname;
    if ( $ref =~ /.*\/([^\/]+)\.[^\.]+$/ ) {
        $refname = $1;
    } else {
        $ref =~ /.*\/([^\/]+)\.?[^\.]*$/;
        $refname = $1;
    }

    # First determine if we are processing bwa_mem or bwa_aln
    if ( $options{'use_mem'} ) {
        print_log_msg( $DEBUG, "Using 'bwa mem'" );

        # If we have 2 paired fastq seqs...
        # Currently not worried about BAM files for 'bwa mem'
        if ($fastq_paired) {
            print_log_msg( $DEBUG, "Detected paired-end fastq" );
            if (   $file_count == 2
                && exists( $type{'fastq_1'} )
                && exists( $type{'fastq_2'} ) )
            {
         # Making the assumption that both FastQ files are in the same directory
                my ( $file_base, $file_dir, $file_ext ) =
                  fileparse( $type{'fastq_1'}, qr/\.[^.]*/ );
                $file_base =~ s/\.fastq$//;
                $out1 = $refname . "_" . $file_base . ".bwa.sam";
				$out1 =~ s/bwa\.prelim\.filtered\.//g;
                my $concat_files = $type{'fastq_1'} . " " . $type{'fastq_2'};
                align_BWA( \%options, "mem", $concat_files, "", $out1, $ref );

            } else {
                print_log_msg( $ERROR,
                    "ERROR : Main :: Irregular number of files $file_count for alignment found in the input directory $options{'input_dir'}"
                );
            }
        } else {
            print_log_msg( $DEBUG, "Detected a single-end fastq sequence" );
            if ( $file_count == 1 && exists( $type{'fastq'} ) ) {
                ( $file_base, $file_dir, $file_ext ) =
                  fileparse( $type{'fastq'}, qr/\.[^.]*/ );
                $file_base =~ s/\.fastq$//;
                $out = $refname . "_" . $file_base . ".bwa.sam";
				$out =~ s/bwa\.prelim\.filtered\.//g;
                align_BWA( \%options, "mem", $type{'fastq'}, "", $out, $ref );
            }
        }
    } else {
        print_log_msg( $DEBUG, "Using 'bwa aln'" );

      # If we have either 2 paired fastq seqs or a single paired-end BAM file...
        if ( $fastq_paired || $options{'bam_paired'} ) {
            print_log_msg( $DEBUG,
                "Detected paired-end fastq or BAM sequences" );
            if ( $file_count == 1 && exists( $type{'bam'} ) ) {
                ( $file_base, $file_dir, $file_ext ) =
                  fileparse( $type{'bam'}, qr/\.[^.]*/ );
                $out1 = $refname . "_" . $file_base . ".aln1.sai";
				$out1 =~ s/bwa\.prelim\.filtered\.//g;
                align_BWA( \%options, "aln", $type{'bam'}, "-b1", $out1, $ref );

                $out2 = $refname . "_" . $file_base . ".aln2.sai";
				$out2 =~ s/bwa\.prelim\.filtered\.//g;
                align_BWA( \%options, "aln", $type{'bam'}, "-b2", $out2, $ref );

                $in =
                    $options{'output_dir'} . "/"
                  . $out1 . " "
                  . $options{'output_dir'} . "/"
                  . $out2 . " "
                  . $type{'bam'} . " "
                  . $type{'bam'};
                $out = $refname . "_" . $file_base . ".bwa.sam";
				$out =~ s/bwa\.prelim\.filtered\.//g;

            } elsif ( $file_count == 2
                && exists( $type{'fastq_1'} )
                && exists( $type{'fastq_2'} ) )
            {
                ( $file_base, $file_dir, $file_ext ) =
                  fileparse( $type{'fastq_1'}, qr/\.[^.]*/ );
                $file_base =~ s/\.fastq$//;
                $out1 = $refname . "_" . $file_base . ".aln1.sai";
				$out1 =~ s/bwa\.prelim\.filtered\.//g;
                align_BWA( \%options, "aln", $type{'fastq_1'}, "", $out1,
                    $ref );

                ( $file_base, $file_dir, $file_ext ) =
                  fileparse( $type{'fastq_2'}, qr/\.[^.]*/ );
                $file_base =~ s/\.fastq$//;
                $out2 = $refname . "_" . $file_base . ".aln2.sai";
				$out2 =~ s/bwa\.prelim\.filtered\.//g;
                align_BWA( \%options, "aln", $type{'fastq_2'}, "", $out2,
                    $ref );

                $file_base =~ s/_2//g;
                $in =
                    $options{'output_dir'} . "/"
                  . $out1 . " "
                  . $options{'output_dir'} . "/"
                  . $out2 . " "
                  . $type{'fastq_1'} . " "
                  . $type{'fastq_2'};
                $out = $refname . "_" . $file_base . ".bwa.sam";
				$out =~ s/bwa\.prelim\.filtered\.//g;

            } else {
                print_log_msg( $ERROR,
                    "ERROR : Main :: Irregular number of files $file_count for alignment found in the input directory $options{'input_dir'}"
                );
            }
            align_BWA( \%options, "sampe", $in, "", $out, $ref );
        } else {
            print_log_msg( $DEBUG,
                "Detected a single-end fastq or BAM sequence" );
            if ( $file_count == 1 && exists( $type{'fastq'} ) ) {
                ( $file_base, $file_dir, $file_ext ) =
                  fileparse( $type{'fastq'}, qr/\.[^.]*/ );
                $file_base =~ s/\.fastq$//;
                $out = $refname . "_" . $file_base . ".aln.sai";
				$out =~ s/bwa\.prelim\.filtered\.//g;
                align_BWA( \%options, "aln", $type{'fastq'}, "", $out, $ref );

                $in =
                  $options{'output_dir'} . "/" . $out . " " . $type{'fastq'};
                $out = $refname . "_" . $file_base . ".bwa.sam";
				$out =~ s/bwa\.prelim\.filtered\.//g;

            } elsif ( $file_count == 1 && exists( $type{'bam'} ) ) {
                ( $file_base, $file_dir, $file_ext ) =
                  fileparse( $type{'bam'}, qr/\.[^.]*/ );
                $out = $refname . "_" . $file_base . ".aln.sai";
				$out =~ s/bwa\.prelim\.filtered\.//g;
                align_BWA( \%options, "aln", $type{'bam'}, "-b0", $out, $ref );

                $in  = $options{'output_dir'} . "/" . $out . " " . $type{'bam'};
                $out = $refname . "_" . $file_base . ".bwa.sam";
				$out =~ s/bwa\.prelim\.filtered\.//g;
            }
            align_BWA( \%options, "samse", $in, "", $out, $ref );
        }
    }

	# Determine if the newly created SAM file is truncated.  Fail if it is
	# Note:  $out is just the basename.  Output_dir is added to command in subroutine
	is_sam_truncated(\%options, $out);
    # Convert SAM file into a BAM file and remove SAM if specified to
    sam_to_bam( \%options, $out );

}
close OUT_LIST;

if ( $options{'keep_sam'} == 0 ) {
    $cmd = "rm "
      . $options{'output_dir'}
      . "/*.sam "
      . $options{'output_dir'}
      . "/*.sai";
    $exit_code = system($cmd);
    if ( $exit_code != 0 ) {
        print_log_msg( $ERROR,
            "ERROR : Main :: Cleanup of SAM and alignment files failed from $options{'output_dir'} with error. Check the stderr"
        );
    } else {
        print_log_msg( $DEBUG,
            "INFO : Main :: Cleanup of SAM and alignment files succesfully completed in $options{'output_dir'}"
        );
    }
}

###############
# SUBROUTINES #
###############

# Convert a SAM-formatted file into a BAM-formatted one
sub sam_to_bam {
    my ( $cmd_line_args, $file ) = @_;
    my ( $cmd, $opts, $file_base, $file_dir, $file_ext, $out );
    my $exit_code;

    my $sub_name = ( caller(0) )[3];
    if ( exists( $cmd_line_args->{'samtools_params'} ) ) {
        $opts = $cmd_line_args->{'samtools_params'};
    } else {
		$opts = '';
	}
    ( $file_base, $file_dir, $file_ext ) = fileparse( $file, qr/\.[^.]*/ );
    $out = $cmd_line_args->{'output_dir'} . "/" . $file_base . ".bam";

    $cmd =
        $cmd_line_args->{'samtools_path'}
      . " view "
      . $opts
      . " -bhS -o "
      . $out . " "
      . $cmd_line_args->{'output_dir'} . "/"
      . $file;
    print_log_msg( $DEBUG,
        "INFO : $sub_name :: Start converting SAM file $cmd_line_args->{'output_dir'}/$file to BAM format $out.\nINFO : $sub_name :: Command : $cmd"
    );
    $exit_code = system($cmd);
    if ( $exit_code != 0 ) {
        print_log_msg( $ERROR,
            "ERROR : $sub_name :: SAM to BAM conversion failed for $cmd_line_args->{'output_dir'}/$file with error. Check the stderr"
        );
    } else {
        print_log_msg( $DEBUG,
            "INFO : $sub_name :: SAM to BAM conversion succesfully completed in $cmd_line_args->{'output_dir'}"
        );
    }

	# Write the new BAM file to our list
	print OUT_LIST $out . "\n";
}

sub is_sam_truncated {
	my $cmd_line_args = shift;
    my $sam = shift;
    my $sub_name = ( caller(0) )[3];
	$cmd = $cmd_line_args->{'samtools_path'} . " view " . $cmd_line_args->{'output_dir'} . "/" . $out;
	print_log_msg($DEBUG, "INFO : $sub_name :: Checking to see if SAM file is truncated\nINFO : $sub_name :: Command : $cmd");
    $exit_code = system($cmd);
    if ( $exit_code != 0 ) {
        print_log_msg( $ERROR,
            "ERROR : $sub_name :: SAM file appears to be truncated.  Try to re-run with a higher memory requirement"
        );
    } else {
        print_log_msg( $DEBUG,
            "INFO : $sub_name :: SAM file seems fine."
        );
    }
	return;
}

sub align_BWA {
    my ( $cmd_line_args, $algo, $files, $opts, $outfile, $ref ) = @_;
    my ( $cmd, $param );
    my $exit_code;
    my %params_hash = (
        'misMsc'  => 'M',
        'maxGapO' => 'o',
        'maxGapE' => 'e',
        'gapOsc'  => 'O',
        'gapEsc'  => 'E',
        'nThrds'  => 't'
    );

    my $sub_name = ( caller(0) )[3];

    if ( $algo eq "aln" ) {
        foreach $param ( keys %params_hash ) {
            if ( exists( $cmd_line_args->{$param} ) ) {
                $opts .=
                  " -" . $params_hash{$param} . " " . $cmd_line_args->{$param};
            }
        }
    } elsif ( $algo eq "mem" ) {
        foreach $param ( keys %params_hash ) {
            if ( exists( $cmd_line_args->{$param} ) ) {
                $opts .=
                  " -" . $params_hash{$param} . " " . $cmd_line_args->{$param};
            }
        }
    } elsif ( $algo eq "sampe" || $algo eq "samse" ) {
        if ( exists( $cmd_line_args->{'maxOcc'} ) ) {
            $opts .= "-n " . $cmd_line_args->{'maxOcc'};
        }
    } else {
        print_log_msg( $ERROR,
            "ERROR : $sub_name :: $algo is not supported by this version of BWA component"
        );
    }

    if ( exists( $cmd_line_args->{'bwa_params'} ) ) {
        $opts .= " " . $cmd_line_args->{'bwa_params'};
    }

	# Same files passed from lgt_bwa component to another lgt_bwa component may have .bwa.prelim.filtered.bam endings.  This doesn't break the script by having it but I would like to keep the name decently short - SAdkins 4/18/16
	$outfile =~ s/bwa\.prelim\.filtered\.//g;

    $cmd =
        $cmd_line_args->{'bwa_path'} . " "
      . $algo . " "
      . $opts . " "
      . $ref . " "
      . $files . " > "
      . $cmd_line_args->{'output_dir'} . "/"
      . $outfile;
    print_log_msg( $DEBUG,
        "INFO : $sub_name :: Start aligning $files to $ref.\nINFO : $sub_name :: Command : $cmd"
    );
    $exit_code = system($cmd);
    if ( $exit_code != 0 ) {
        print_log_msg( $ERROR,
            "ERROR : $sub_name :: $files alignment failed with error. Check the stderr"
        );
    } else {
        print_log_msg( $DEBUG,
            "INFO : $sub_name :: $files alignment to $ref succesfully completed in $cmd_line_args->{'output_dir'}"
        );
    }
}

# Determine format of input file or files from input directory to be aligned
sub determine_format {
    my ( $opts, $file_type ) = @_;
    my $file = $opts->{'input_file'};

    # Going to check and see if we have a inputted file first
    if ($file) {
        my ( $base, $dir_path, $ext ) = fileparse( $file, qr/\.[^.]*/ );

        # Depending on the extension, we handle it differently
        if ( $ext =~ /list/ ) {
            my @list_files = `cat $file`;

# SAdkins - 7/7/16  the filter_dups_lc_seqs component will now output multiple BAM files since it is iterative. Because of this I need to merge the BAM files from the list into a single BAM file.
            if ( scalar @list_files > 1 && $list_files[0] =~ /\.bam$/ ) {
            	print_log_msg($DEBUG, "Found multiple BAM files in list.  Merging into single BAM input");
                my $merged_bam = merge_bam_files($file, $opts);
                @list_files = ();
                push @list_files, $merged_bam;
            }

# Currently the script cannot process multiple samples.  Ideally the list file should just contain 1 sample, which is either 1 single-end fastq, 1 BAM, or 2 paired-end fastq files.  Also the .blank file originating from the sra2fastq component in Ergatis can exist in a list file
            foreach my $f (@list_files) {
				chomp $f;

                my ( $base, $dir_path, $ext ) = fileparse( $f, qr/\.[^.]*/ );
                if ( $ext =~ /fastq$/ ) {
                    if ( $ext =~ /_1\./ ) {
                        set_query_to_file_type( $file_type, 'fastq_1', $f );
                        $fastq_paired = 1;
                    } elsif ( $ext =~ /_2\./ ) {
                        set_query_to_file_type( $file_type, 'fastq_2', $f );
                    } else {
                        set_query_to_file_type( $file_type, 'fastq', $f );
                        $fastq_paired = 0;
                    }
                } else {

          # For the BAM, or .blank extensions, we can just recycle current code.
                    my $temp_opts = {};
                    $temp_opts->{'input_file'} = $f;
                    determine_format( $temp_opts, $file_type );
                }
            }
        } elsif ( $ext =~ /bam/ ) {
            set_query_to_file_type( $file_type, 'bam', $file );
        } elsif ( $ext =~ /fastq/ ) {    # Single-end FASTQ file
            set_query_to_file_type( $file_type, 'fastq', $file );
            $fastq_paired = 0;
        } elsif ( $ext =~ /blank/ ) {

# My way of grouping paired FASTQ files is to use the basename in a .blank file
# .blank is just to indicate the 2 paired-end fastq files exist in that directory
            grab_files_from_dir( $dir_path, $file_type );
        }
    } else {

       # If an inputted file wasn't passed, then it has to be an input directory
        my $dir = $opts->{'input_dir'};
        grab_files_from_dir( $dir, $file_type );
    }
}

# Read through a directory to grab relevant fastq or bam files
sub grab_files_from_dir {
    my ( $dir, $file_type ) = @_;
    my $file;
    my $sub_name = ( caller(0) )[3];
    opendir( DIR, $dir )
      or print_log_msg( $ERROR,
        "ERROR : $sub_name :: Could not open directory $dir for reading.\nReason : $!"
      );

      # TODO: write merged bam files subroutine for directory globbing
    while ( $file = readdir(DIR) ) {
        my $path = $dir . "/" . $file;
        if ( $file =~ /fastq$/ ) {
            if ( $file =~ /_1\./ ) {
                set_query_to_file_type( $file_type, 'fastq_1', $path );
                $fastq_paired = 1;
            } elsif ( $file =~ /_2\./ ) {
                set_query_to_file_type( $file_type, 'fastq_2', $path );
            } else {
                set_query_to_file_type( $file_type, 'fastq', $path );
                $fastq_paired = 0;
            }
        } elsif ( $file =~ /bam$/ ) {
            set_query_to_file_type( $file_type, 'bam', $path );
        }
    }
}

# This function accepts the input BAM list file and uses it as an argument for
#   Samtools merge to merge all BAM files in the list into a single BAM file
sub merge_bam_files {
    my $sub_name = ( caller(0) )[3];
    my $list_file = shift;
	my $cmd_line_args = shift;
    my $merged_bam = $cmd_line_args->{'output_dir'} . "/" . $$ ."_merged.bam";
    my $cmd = $cmd_line_args->{'samtools_path'} . " merge -b " . $list_file . " " . $merged_bam;
	print_log_msg($DEBUG, "INFO : $sub_name :: Command : $cmd");
    $exit_code = system($cmd);
    if ( $exit_code != 0 ) {
        print_log_msg( $ERROR,
            "ERROR : $sub_name :: Merging BAM files failed for $merged_bam with error. Check the stderr"
        );
    } else {
        print_log_msg( $DEBUG,
            "INFO : $sub_name :: BAM file merging succesfully completed for $merged_bam"
        );
    }

    return $merged_bam;
}

# This function checks to see if a file has already been associated with the given type.
# If not, then assigns that file path to that type in the hash.
# If so, errors out, since only one query sample is allowed.
sub set_query_to_file_type {
    my ( $type_h, $type, $path ) = @_;

    print_log_msg( $ERROR,
        "ERROR : File path of the type $type has already been added.  Please just provide one query sample input for alignment."
    ) if ( exists $type_h->{$type} );
    $type_h->{$type} = $path;
}

# This function ensures that only one sample query was provided.
sub check_for_single_sample {
    my $type = shift;
    if ( exists $type->{'bam'} ) {
        if (   exists $type->{'fastq'}
            || exists $type->{'fastq_1'}
            || exists $type->{'fastq_2'} )
        {
            print_log_msg( $ERROR,
                "ERROR : Detected both fastq and BAM files as query input. Please just provide only one sample query input for alignment."
            );
        }
    } elsif ( exists $type->{'fastq'} ) {
        if ( exists $type->{'fastq_1'} || exists $type->{'fastq_2'} ) {
            print_log_msg( $ERROR,
                "ERROR : Detected both single and paired-end FASTQ files as query input. Please just provide only one sample query input for alignment."
            );
        }
    }
}

# Description   : Used to check the correctness of the command line arguments passed to the script. The script exits if required arguments are missing.
# Parameters    : NA
# Returns       : NA
# Modifications :

sub check_parameters {
    my ($opts) = @_;
    my $sub_name = ( caller(0) )[3];
    if ( exists( $opts->{'log'} ) ) {
        open( $log_fh, "> $opts->{'log'}" )
          or die "Could not open $opts->{'log'} file for writing.Reason : $!\n";
    }
    my @required = qw( output_dir bwa_path);
    foreach my $opt (@required) {
        if ( !defined( $opts->{$opt} ) ) {
            print_log_msg( $ERROR,
                "ERROR : $sub_name :: Required option $opt not passed" );
        }
    }

    if ( $opts->{'reference'} ) {
        @ref_files = split( /,/, $opts->{'reference'} );
    } elsif ( $opts->{'reference_list'} ) {
        @ref_files = `cat $opts->{'reference_list'}`;
    } else {
        print_log_msg( $ERROR,
            "Either --reference or --reference_list must be passed" );
    }

    print_log_msg( $ERROR,
        "ERROR : $sub_name :: Either --input_file or --input_dir are required" )
      if ( !(defined $opts->{'input_file'} || defined $opts->{'input_dir'}) );

	# Check for empty file for input file
	if (defined $opts->{'input_file'} && !-s $opts->{'input_file'}) {
		print_log_msg( $ERROR, "ERROR : $sub_name :: Input file or list " . $opts->{'input_file'} . " appears to be empty" )
	}

    # Delete SAM files once they are converted to BAM by default
    if ( !defined( $opts->{'keep_sam'} ) ) {
        $opts->{'keep_sam'} = 0;
    }
}

####################################################################################################################################################

# Description   : Used to handle logging of messages(errors and warnings) during the execution of the script
# Parameters    : level = can be ERROR, WARNING or INFO
#		  msg   = msg to be printed in the log file or to STDERR
# Returns       : NA
# Modifications :

sub print_log_msg {
    my ( $level, $msg ) = @_;
    if ( $level <= $DEBUG ) {
        print STDERR "$msg\n";
        die "" if ( $level == $ERROR );
    }
    print $log_fh "$msg\n" if ( defined($log_fh) );
}

__END__

#####################
# POD DOCUMENTATION #
#####################

=head1 NAME

# Name of the script and a 1 line desc

=head1 SYNOPSIS

# USAGE :

	parameters in [] are optional

=head1 OPTIONS

B<--reference,-r>
	Name of the fasta reference to align against.  Multiple files can be specified if seperated by commas

B<--reference_list, -R>
	List of fasta references to align against.  The list is looped over and each reference aligned to individually.

B<--input_file, -i>
	Fastq or BAM input file.  A list file pointing to input is also accepted.

B<--input_dir,-I>
	Input directory where either fastq files or BAM input files can be found

B<--output_dir, -o>
	Directory to store the output contents

B<--misMsc>
	(Optional) Mismatch penalty. BWA will not search for suboptimal hits with a score lower than (bestScore-misMsc).  Default is 3

B<--maxGapO>
    (Optional)  Maximum number of gap opens.  Default is 1

B<--maxGapE>
	(Optional) Maximum number of gap extensions, -1 for k-difference mode (disallowing long gaps).  Default is -1

B<--maxOcc>
	(Optional) Maximum edit distance if the value is INT, or the fraction of missing alignments given 2% uniform base error rate if FLOAT. In the latter case, the maximum edit distance is automatically chosen for different read lengths. Default is 0.04

B<--gapOsc>
	(Optional) Gap open penalty. Default is 11

B<--gapEsc>
	(Optional) Gap extension penalty. Default is 4

B<--nThrds, -t>
	(Optional) Number of threads used in multi-threading mode.  Not necessary if calling via Ergatis. Default is 1

B<--bwa_params, -B>
	Extra parameters for the "bwa" program

B<--samtools_params, -S>
	Extra parameters for the "samtools" program

B<--bwa_path, -b>
	Path to the "bwa" executable

B<--samtools_path, -s>
	Path to the "samtools" executable

B<--keep_sam, -k>
	If set to 1, SAM files are kept after converting to BWA. Set to 0 and they are removed from the output directory

B<--use_mem, -m>
	Set to 1 to indicate that the "bwa mem" algorithm should be used instead of "bwa aln"

B<--bam_paired, -p>
	If set to 1, this indicates the input BAM is from paired-end data

B<--log, -l>
	Path to writable log file

B<--help, -h>
	Prints this help doc

=head1 DESCRIPTION



=head1 INPUT



=head1 OUTPUT



=head1 AUTHOR

	Shaun Adkins
	Bioinformatics Software Engineer
	Institute for Genome Sciences
	School of Medicine
	University of Maryland
	Baltimore, MD - 21201
	sadkins@som.umaryland.edu

==cut
