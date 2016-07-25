#!/usr/bin/perl -w

eval 'exec /usr/bin/perl -w -S $0 ${1+"$@"}'
    if 0; # not running under some shell

#########################################################################################
#											#
# Name	      : lgt_filter_reads.pl							#
# Version     : 1.0									#
# Project     : LGT Seek Pipeline							#
# Description : Script to filter out non LGT aligned reads 				#
# Author      : Sonia Agrawal								#
# Date        : September 25, 2015							#
#											#
#########################################################################################

use strict;

# Extended processing of command line options
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);

# Prints usage message from embedded pod documentation
use Pod::Usage;

# Include packages here
use File::Basename;
use LGT::Common;

#############
# CONSTANTS #
#############

###########
# GLOBALS #
###########
my %hCmdLineArgs = ();
my @input_files  = [];

# Log file handle;
my $fhLog;
my ( $ERROR, $WARN, $DEBUG ) = ( 1, 2, 3 );
my %hType = ();
my ( $nFileCnt, $nExitCode );
my ( $sIn, $sOut, $sOut1, $sOut2, $sFbase, $sFdir, $sFext, $sCmd, $sOutFile );
my $iRead = 0;
my ( $sQnameR1, $iFlagR1, $sCigarR1, $sQnameR2, $iFlagR2, $sCigarR2 );
my %hStatR1     = ();
my %hStatR2     = ();
my $iSingletons = 0;
my $num_null    = 0;
my $SC          = 0;

################
# MAIN PROGRAM #
################
GetOptions(
    \%hCmdLineArgs,      'input_file|i=s',
    'input_list|I=s',    'output_dir|o=s',
    'samtools_path|s=s', 'samtools_params|S=s',
    'softclip_min|m=i',  'keep_mm_reads|mm=i',
    'log|l=s',           'help|h'
) or pod2usage();

pod2usage( { -exitval => 0, -verbose => 2, -output => \*STDERR } )
  if ( $hCmdLineArgs{'help'} );

checkCmdLineArgs( \%hCmdLineArgs );

foreach my $i_file (@input_files) {
    chomp $i_file;

    # MM, MU, and UU, keep count of mapped and unmapped mate pairs per BAM input
    my $MM = 0;
    my $MU = 0;    # Also keeps track of UM
    my $UU = 0;
    ( $sFbase, $sFdir, $sFext ) = fileparse( $i_file, qr/\.[^.]*/ );
    my $sSortedFile = $sFbase . ".name.sorted.bam";
    SortBam( \%hCmdLineArgs, $sSortedFile, $i_file );

    # View the sorted BAM file, and only output the SAM headers
    $sCmd =
        $hCmdLineArgs{'samtools_path'}
      . " view -H "
      . $hCmdLineArgs{'output_dir'} . "/"
      . $sSortedFile;
    my $sHeader = `$sCmd`;
    if ( $? != 0 ) {
        printLogMsg( $ERROR,
            "ERROR : main :: Retrieving of the header from the sorted BAM $hCmdLineArgs{'output_dir'}"
              . "/"
              . $sSortedFile
              . " failed. Check the stderr" );
    }

    $sOutFile =
      $hCmdLineArgs{'output_dir'} . "/" . $sFbase . ".prelim.filtered.bam";

    # print SAM headers into a new BAM file
    open( my $fhFW, "| $hCmdLineArgs{'samtools_path'} view -S - -bo $sOutFile" )
      or printLogMsg(
        $ERROR,
        "ERROR : main :: Could not open BAM file $sOutFile for writing.\nReason : $!"
      );
    print $fhFW $sHeader;

    # Open the sorted BAM file for reading
    open(
        my $fhFR,
        "$hCmdLineArgs{'samtools_path'} view "
          . $hCmdLineArgs{'output_dir'} . "/"
          . $sSortedFile . " |"
      )
      or printLogMsg(
        $ERROR,
        "ERROR : main :: Could not open BAM file "
          . $hCmdLineArgs{'output_dir'} . "/"
          . $sSortedFile
          . " for reading.\nReason : $!"
      );
    my $sRead1;

    while ( my $sLine = <$fhFR> ) {
        chomp($sLine);

        # Keeping track if current line is first or second mate pair
        if ( $iRead == 0 ) {
            $sRead1 = $sLine;
            $iRead  = 1;
            next;
        }
        $iRead = 0;
        my $print  = 0;
        my $sRead2 = $sLine;
        ( $sQnameR1, $iFlagR1, $sCigarR1 ) = ( split /\t/, $sRead1 )[ 0, 1, 5 ];
        ( $sQnameR2, $iFlagR2, $sCigarR2 ) = ( split /\t/, $sRead2 )[ 0, 1, 5 ];

        # Ensure the two reads actually have valid names
        if ( $sQnameR1 =~ /null/ ) {
            $num_null++;

            # Treat second read as the new first, and grab a new second.
            $sRead1 = $sRead2;
            $iRead  = 1;
            next;
        }
        if ( $sQnameR2 =~ /null/ ) {
            $num_null++;
            $iRead = 1;
            next;
        }

        # Ensure the two reads are actually mates
        if ( $sQnameR1 ne $sQnameR2 ) {
            $iSingletons++;

     # If not mates, then treat second read as the first, and grab a new second.
            $sRead1 = $sRead2;
            $iRead  = 1;
            next;
        }

        my $stat_r1 = parse_flag($iFlagR1);
        my $stat_r2 = parse_flag($iFlagR2);

# Want to keep UU, MU, and UM reads.  Keep MM if specified (reference genome is donor for example)
        if ( !$stat_r1->{'qunmapped'} && !$stat_r2->{'qunmapped'} ) {
            $print = 1 if $hCmdLineArgs{'keep_mm_reads'};
            $MM += 2;
        }
        if ( $stat_r1->{'qunmapped'} || $stat_r2->{'qunmapped'} ) {

            # Print unmapped reads to filtered SAM file
            $print = 1;

            if ( !$stat_r1->{'qunmapped'} && $stat_r2->{'qunmapped'} ) {
                $MU += 2;
            }
            if ( $stat_r1->{'qunmapped'} && !$stat_r2->{'qunmapped'} ) {
                $MU += 2;
            }
            if ( $stat_r1->{'qunmapped'} && $stat_r2->{'qunmapped'} ) {
                $UU += 2;
            }

        }

        ## FILTER FOR soft clipped reads
        if ( $hCmdLineArgs{'softclip_min'} ) {
            map {
                if ( $_ =~ /(\d+)S/ && $1 >= $hCmdLineArgs{'softclip_min'} ) {
                    $print = 1;
                    $SC += 2;
                }    ## 10.21.14 Testing more liberal softclip parsing
            } ( $sCigarR1, $sCigarR2 );
        }

        if ($print) {
            print $fhFW $sRead1 . "\n";
            print $fhFW $sRead2 . "\n";
        }
    }
    close $fhFW;
    printLogMsg( $DEBUG,
        "INFO : Main :: BAM Filtering results for $sFbase:\n\tBad Singletons : $iSingletons\n\tNull : $num_null\n\tSoft-clipped Reads : $SC\n\tMapped-Unmapped : $MU\n\tMapped-Mapped : $MM\n\tUnmapped-Unmapped : $UU\n"
    );
}

###############
# SUBROUTINES #
###############

# Sort BAM file by read names
sub SortBam {
    my ( $phCmdLineArgs, $sFile, $input_file ) = @_;
    my ( $sCmd, $sOut );
    my $nExitCode;
    my $sOptions = "";

    my $sSubName = ( caller(0) )[3];
    if ( exists( $phCmdLineArgs->{'samtools_params'} ) ) {
        $sOptions = $phCmdLineArgs->{'samtools_params'};
    }
    $sOut = $phCmdLineArgs->{'output_dir'} . "/" . $sFile;

$sCmd = $phCmdLineArgs->{'samtools_path'}." sort ".$sOptions." -O bam -n -o ".$sOut." -T /tmp/".$sFile." " . $input_file;
    printLogMsg( $DEBUG,
        "INFO : $sSubName :: Start sorting file $input_file to produce sorted BAM $sOut.\nINFO : $sSubName :: Command : $sCmd"
    );
    $nExitCode = system($sCmd);
    if ( $nExitCode != 0 ) {
        printLogMsg( $ERROR,
            "ERROR : $sSubName :: Sorting of BAM failed for $input_file with error. Check the stderr"
        );
    } else {
        printLogMsg( $DEBUG,
            "INFO : $sSubName :: Sorting of BAM succesfully completed in $phCmdLineArgs->{'output_dir'}/"
              . $sFile );
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
    my @aRequired = qw( output_dir samtools_path);
    foreach my $sOption (@aRequired) {
        if ( !defined( $phCmdLineArgs->{$sOption} ) ) {
            printLogMsg( $ERROR,
                "ERROR : $sSubName :: Required option $sOption not passed" );
        }
    }

    if ( $phCmdLineArgs->{'input_file'} ) {
        push @input_files, $phCmdLineArgs->{'input_file'};
    } elsif ( $phCmdLineArgs->{'input_list'} ) {
        @input_files = `cat $phCmdLineArgs->{'input_list'}`;
    } else {
        printLogMsg( $ERROR,
            "ERROR : Either --input_file or --input_file_list must be passed" );
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
