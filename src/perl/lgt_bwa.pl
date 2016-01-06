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
my %hCmdLineArgs = ();

# Log file handle;
my $fhLog;
my ( $ERROR, $WARN, $DEBUG ) = ( 1, 2, 3 );
my %hType = ();
my ( $nFileCnt, $nExitCode );
my ( $sIn, $sOut, $sOut1, $sOut2, $sFbase, $sFdir, $sFext, $sCmd );
my $fastq_paired = 0;    # Determines if fastq seqs are paired-end or not

################
# MAIN PROGRAM #
################
GetOptions(
    \%hCmdLineArgs,        'reference|r=s',
    'bam_paired|p=i',      'input_dir|I=s',
    'input_file|i=s',      'output_dir|d=s',
    'bwa_path|b=s',        'samtools_path|s=s',
    'misMsc|M=i',          'maxGapO|o=i',
    'maxGapE|e=i',         'gapOsc|O=i',
    'gapEsc|E=i',          'nThrds|t=i',
    'maxOcc|n=i',          'bwa_params|B=s',
    'samtools_params|S=s', 'keep_sam|k=i',
    'log|l=s',             'help|h'
) or pod2usage();

pod2usage( { -exitval => 0, -verbose => 2, -output => \*STDERR } )
  if ( $hCmdLineArgs{'help'} );

checkCmdLineArgs( \%hCmdLineArgs );

#Go though the directory and determine the type of files that are to be aligned
DetermineFormat( \%hCmdLineArgs, \%hType );

$nFileCnt = keys %hType;

# If we have either 2 paired fastq seqs or a single paired-end BAM file...
if ( $fastq_paired || $hCmdLineArgs{'bam_paired'} ) {
    printLogMsg( $DEBUG, "Detected paired-end fastq or BAM sequences" );
    if ( $nFileCnt == 1 && exists( $hType{'bam'} ) ) {
        ( $sFbase, $sFdir, $sFext ) = fileparse( $hType{'bam'}, qr/\.[^.]*/ );
        $sOut1 = $sFbase . ".aln1.sai";
        AlignBwa( \%hCmdLineArgs, "aln", $hType{'bam'}, "-b1", $sOut1 );

        $sOut2 = $sFbase . ".aln2.sai";
        AlignBwa( \%hCmdLineArgs, "aln", $hType{'bam'}, "-b2", $sOut2 );

        $sIn =
            $hCmdLineArgs{'output_dir'} . "/"
          . $sOut1 . " "
          . $hCmdLineArgs{'output_dir'} . "/"
          . $sOut2 . " "
          . $hType{'bam'} . " "
          . $hType{'bam'};
        $sOut = $sFbase . ".bwa.sam";

        #		AlignBwa(\%hCmdLineArgs, "sampe", $sIn , "", $sOut);

    } elsif ( $nFileCnt == 2
        && exists( $hType{'fastq_1'} )
        && exists( $hType{'fastq_2'} ) )
    {
        ( $sFbase, $sFdir, $sFext ) =
          fileparse( $hType{'fastq_1'}, qr/\.[^.]*/ );
        $sFbase =~ s/\.fastq$//;
        $sOut1 = $sFbase . ".aln1.sai";
        AlignBwa( \%hCmdLineArgs, "aln", $hType{'fastq_1'}, "", $sOut1 );

        ( $sFbase, $sFdir, $sFext ) =
          fileparse( $hType{'fastq_2'}, qr/\.[^.]*/ );
        $sFbase =~ s/\.fastq$//;
        $sOut2 = $sFbase . ".aln2.sai";
        AlignBwa( \%hCmdLineArgs, "aln", $hType{'fastq_2'}, "", $sOut2 );

        $sFbase =~ s/_2//g;
        $sIn =
            $hCmdLineArgs{'output_dir'} . "/"
          . $sOut1 . " "
          . $hCmdLineArgs{'output_dir'} . "/"
          . $sOut2 . " "
          . $hType{'fastq_1'} . " "
          . $hType{'fastq_2'};
        $sOut = $sFbase . ".bwa.sam";

        #		AlignBwa(\%hCmdLineArgs, "sampe", $sIn , "", $sOut);

    } else {
        printLogMsg( $ERROR,
            "ERROR : Main :: Irregular number of files $nFileCnt for alignment found in the input directory $hCmdLineArgs{'input_dir'}"
        );
    }
    AlignBwa( \%hCmdLineArgs, "sampe", $sIn, "", $sOut );
} else {
    printLogMsg( $DEBUG, "Detected a single-end fastq or BAM sequence" );
    if ( $nFileCnt == 1 && exists( $hType{'fastq'} ) ) {
        ( $sFbase, $sFdir, $sFext ) = fileparse( $hType{'fastq'}, qr/\.[^.]*/ );
        $sFbase =~ s/\.fastq$//;
        $sOut = $sFbase . ".aln.sai";
        AlignBwa( \%hCmdLineArgs, "aln", $hType{'fastq'}, "", $sOut );

        $sIn =
          $hCmdLineArgs{'output_dir'} . "/" . $sOut . " " . $hType{'fastq'};
        $sOut = $sFbase . ".bwa.sam";

        #		AlignBwa(\%hCmdLineArgs, "samse", $sIn , "", $sOut);

    } elsif ( $nFileCnt == 1 && exists( $hType{'bam'} ) ) {
        ( $sFbase, $sFdir, $sFext ) = fileparse( $hType{'bam'}, qr/\.[^.]*/ );
        $sOut = $sFbase . ".aln.sai";
        AlignBwa( \%hCmdLineArgs, "aln", $hType{'bam'}, "-b0", $sOut );

        $sIn  = $hCmdLineArgs{'output_dir'} . "/" . $sOut . " " . $hType{'bam'};
        $sOut = $sFbase . ".bwa.sam";

        #		AlignBwa(\%hCmdLineArgs, "samse", $sIn , "", $sOut);
    }
    AlignBwa( \%hCmdLineArgs, "samse", $sIn, "", $sOut );
}

SamtoBam( \%hCmdLineArgs, $sOut );
if ( $hCmdLineArgs{'keep_sam'} == 0 ) {
    $sCmd = "rm "
      . $hCmdLineArgs{'output_dir'}
      . "/*.sam "
      . $hCmdLineArgs{'output_dir'}
      . "/*.sai";
    $nExitCode = system($sCmd);
    if ( $nExitCode != 0 ) {
        printLogMsg( $ERROR,
            "ERROR : Main :: Cleanup of SAM and alignment files failed from $hCmdLineArgs{'output_dir'} with error. Check the stderr"
        );
    } else {
        printLogMsg( $DEBUG,
            "INFO : Main :: Cleanup of SAM and alignment files succesfully completed in $hCmdLineArgs{'output_dir'}"
        );
    }
}

###############
# SUBROUTINES #
###############

sub SamtoBam {
    my ( $phCmdLineArgs, $sFile ) = @_;
    my ( $sCmd, $sOptions, $sFbase, $sFdir, $sFext, $sOut );
    my $nExitCode;

    my $sSubName = ( caller(0) )[3];
    if ( exists( $phCmdLineArgs->{'samtools_params'} ) ) {
        $sOptions = $phCmdLineArgs->{'samtools_params'};
    }
    ( $sFbase, $sFdir, $sFext ) = fileparse( $sFile, qr/\.[^.]*/ );
    $sOut = $phCmdLineArgs->{'output_dir'} . "/" . $sFbase . ".bam";

    $sCmd =
        $phCmdLineArgs->{'samtools_path'}
      . " view "
      . $sOptions
      . " -bhS -o "
      . $sOut . " "
      . $phCmdLineArgs->{'output_dir'} . "/"
      . $sFile;
    printLogMsg( $DEBUG,
        "INFO : $sSubName :: Start converting SAM file $phCmdLineArgs->{'output_dir'}/$sFile to BAM format $sOut.\nINFO : $sSubName :: Command : $sCmd"
    );
    $nExitCode = system($sCmd);
    if ( $nExitCode != 0 ) {
        printLogMsg( $ERROR,
            "ERROR : $sSubName :: SAM to BAM conversion failed for $phCmdLineArgs->{'output_dir'}/$sFile with error. Check the stderr"
        );
    } else {
        printLogMsg( $DEBUG,
            "INFO : $sSubName :: SAM to BAM conversion succesfully completed in $phCmdLineArgs->{'output_dir'}"
        );
    }
}

sub AlignBwa {
    my ( $phCmdLineArgs, $sAlgo, $sFiles, $sOptions, $sOutFile ) = @_;
    my ( $sCmd, $sParam );
    my $nExitCode;
    my %hParams = (
        'misMsc'  => 'M',
        'maxGapO' => 'o',
        'maxGapE' => 'e',
        'gapOsc'  => 'O',
        'gapEsc'  => 'E',
        'nThrds'  => 't'
    );

    my $sSubName = ( caller(0) )[3];

    if ( $sAlgo eq "aln" ) {
        foreach $sParam ( keys %hParams ) {
            if ( exists( $phCmdLineArgs->{$sParam} ) ) {
                $sOptions .=
                  " -" . $hParams{$sParam} . " " . $phCmdLineArgs->{$sParam};
            }
        }
    } elsif ( $sAlgo eq "sampe" || $sAlgo eq "samse" ) {
        if ( exists( $phCmdLineArgs->{'maxOcc'} ) ) {
            $sOptions .= "-n " . $phCmdLineArgs->{'maxOcc'};
        }
    } else {
        printLogMsg( $ERROR,
            "ERROR : $sSubName :: $sAlgo is not supported by this version of BWA component"
        );
    }

    if ( exists( $phCmdLineArgs->{'bwa_params'} ) ) {
        $sOptions .= " " . $phCmdLineArgs->{'bwa_params'};
    }

	my $ref = $phCmdLineArgs->{'reference'};
	$sCmd =
	    $phCmdLineArgs->{'bwa_path'} . " "
	    . $sAlgo . " "
	    . $sOptions . " "
	    . $ref . " "
	    . $sFiles . " > "
	    . $phCmdLineArgs->{'output_dir'} . "/"
	    . $sOutFile;
	printLogMsg( $DEBUG,
	    "INFO : $sSubName :: Start aligning $sFiles to $ref.\nINFO : $sSubName :: Command : $sCmd"
	);
	$nExitCode = system($sCmd);
	if ( $nExitCode != 0 ) {
	    printLogMsg( $ERROR,
	        "ERROR : $sSubName :: $sFiles alignment failed with error. Check the stderr"
	    );
	} else {
	    printLogMsg( $DEBUG,
	        "INFO : $sSubName :: $sFiles alignment to $ref succesfully completed in $phCmdLineArgs->{'output_dir'}"
	    );
	}
}

# Determine format of input file or files from input directory to be aligned
sub DetermineFormat {
    my ( $phCmdLineArgs, $phType ) = @_;
    my $file = $phCmdLineArgs->{'input_file'};

    # Going to check and see if we have a inputted file first
    if ($file) {
        my ( $base, $dir_path, $ext ) = fileparse( $file, qr/\.[^.]*/ );
        if ( $ext =~ /bam/ ) {
            $phType->{'bam'} = $file;
            return;
        } elsif ( $ext =~ /fastq/ ) {    # Single-end FASTQ file
            $phType->{'fastq'} = $file;
            $fastq_paired = 0;
            return;
        } elsif ( $ext =~ /blank/ )
        { # My way of grouping paired FASTQ files is to use the basename in a .blank file
            GrabFilesFromDir( $dir_path, $phType );
        }
    } else {

       # If an inputted file wasn't passed, then it has to be an input directory
        my $dir = $phCmdLineArgs->{'input_dir'};
        GrabFilesFromDir( $dir, $phType );
    }
}

# Read through a directory to grab relevant fastq or bam files
sub GrabFilesFromDir {
    my ( $dir, $file_type ) = @_;
    my $sFile;
    my $sSubName = ( caller(0) )[3];
    opendir( DIR, $dir )
      or printLogMsg( $ERROR,
        "ERROR : $sSubName :: Could not open directory $dir for reading.\nReason : $!"
      );
    while ( $sFile = readdir(DIR) ) {
        my $sPath = $dir . "/" . $sFile;
        if ( $sFile =~ /fastq$/ ) {
            if ( $sFile =~ /_1\./ ) {
                $file_type->{'fastq_1'} = $sPath;
                $fastq_paired = 1;
            } elsif ( $sFile =~ /_2\./ ) {
                $file_type->{'fastq_2'} = $sPath;
            } else {
                $file_type->{'fastq'} = $sPath;
                $fastq_paired = 0;
            }
        } elsif ( $sFile =~ /bam$/ ) {
            $file_type->{'bam'} = $sPath;
        }
    }
}

# Description   : Used to check the correctness of the command line arguments passed to the script. The script exits if required arguments are missing.
# Parameters    : NA
# Returns       : NA
# Modifications :

sub checkCmdLineArgs {
    my ($phCmdLineArgs) = @_;
    my $sSubName = ( caller(0) )[3];
    if ( exists( $phCmdLineArgs->{'log'} ) ) {
        open( $fhLog, "> $phCmdLineArgs->{'log'}" )
          or die
          "Could not open $phCmdLineArgs->{'log'} file for writing.Reason : $!\n";
    }
    my @aRequired = qw(reference output_dir bwa_path);
    foreach my $sOption (@aRequired) {
        if ( !defined( $phCmdLineArgs->{$sOption} ) ) {
            printLogMsg( $ERROR,
                "ERROR : $sSubName :: Required option $sOption not passed" );
        }
    }

    printLogMsg( $ERROR,
        "ERROR : $sSubName :: Either --input_file or --input_dir are required" )
      if ( !defined $phCmdLineArgs->{'input_dir'}
        && !defined $phCmdLineArgs->{'input_file'} );

    # Delete SAM files once they are converted to BAM by default
    if ( !defined( $phCmdLineArgs->{'keep_sam'} ) ) {
        $phCmdLineArgs->{'keep_sam'} = 0;
    }
}

####################################################################################################################################################

# Description   : Used to handle logging of messages(errors and warnings) during the execution of the script
# Parameters    : level = can be ERROR, WARNING or INFO
#		  msg   = msg to be printed in the log file or to STDERR
# Returns       : NA
# Modifications :

sub printLogMsg {
    my ( $nLevel, $sMsg ) = @_;
    if ( $nLevel <= $DEBUG ) {
        print STDERR "$sMsg\n";
        die "" if ( $nLevel == $ERROR );
    }
    print $fhLog "$sMsg\n" if ( defined($fhLog) );
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



=head1 DESCRIPTION



=head1 INPUT



=head1 OUTPUT



=head1 AUTHOR

	Sonia Agrawal
	Bioinformatics Software Engineer II
	Institute for Genome Sciences
	School of Medicine
	University of Maryland
	Baltimore, MD - 21201
	sagrawal@som.umaryland.edu

==cut
