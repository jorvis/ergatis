#!/usr/local/bin/perl -w

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

#############
# CONSTANTS #
#############


###########
# GLOBALS #
###########
my %hCmdLineArgs = ();
# Log file handle;
my $fhLog;
my ($ERROR, $WARN, $DEBUG) = (1, 2, 3);
my %hType = ();
my ($nFileCnt, $nExitCode);
my ($sIn, $sOut, $sOut1, $sOut2, $sFbase, $sFdir, $sFext, $sCmd, $sOutFile);
my $iRead = 0;
my ($sQnameR1, $iFlagR1, $sCigarR1, $sQnameR2, $iFlagR2, $sCigarR2);
my %hStatR1 = ();
my %hStatR2 = ();

################
# MAIN PROGRAM #
################
GetOptions(\%hCmdLineArgs,
	   'input_file|i=s',
	   'output_dir|o=s',
	   'samtools_path|s=s',
	   'samtools_params|S=s',
	   'log|l=s', 
	   'help|h'
	  ) or pod2usage();

pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDERR} ) if ($hCmdLineArgs{'help'});

checkCmdLineArgs(\%hCmdLineArgs);

($sFbase,$sFdir,$sFext) = fileparse($hCmdLineArgs{'input_file'}, qr/\.[^.]*/);
$sSortedFile = $sFbase."name.sorted.bam";

SortBam(\%hCmdLineArgs, $sSortedFile);

$sCmd = $hCmdLineArgs{'samtools_path'}." view -H ".$hCmdLineArgs{'output_dir'}."/".$sSortedFile;
$sHeader = `$sCmd`;
if($? != 0) {
	printLogMsg($ERROR, "ERROR : main :: Retrieving of the header from the sorted BAM $hCmdLineArgs{'output_dir'}"."/".$sSortedFile." failed. Check the stderr");	
}

$sOutFile = $hCmdLineArgs{'output_dir'}."/".$sFbase.".prelim.filtered.bam";

open($fhFW, "| $hCmdLineArgs{'samtools_path'} view -S - -bo $sOutFile") or printLogMsg($ERROR, "ERROR : main :: Could not open BAM file $sOutFile for writing.\nReason : $!"); 

print $fhFW $sHeader;

open($fhFR, "$hCmdLineArgs{'samtools_path'} view ".$hCmdLineArgs{'output_dir'}."/".$sSortedFile." |") or printLogMsg($ERROR, "ERROR : main :: Could not open BAM file ".$hCmdLineArgs{'output_dir'}."/".$sSortedFile." for reading.\nReason : $!");

while($sLine = <$fhFR>) {
	chomp($sLine);
	%hRStat1 = ();
	%hRStat2 = ();
	if($iRead == 0) {
		$sRead1 = $sLine;
		$iRead = 1;
		next;
	}
	$iRead = 0;
	$sRead2 = $sLine;
	($sQnameR1, $iFlagR1, $sCigarR1) = (split /\t/, $sRead1)[ 0, 1, 5];
	($sQnameR2, $iFlagR2, $sCigarR2) = (split /\t/, $sRead2)[ 0, 1, 5];
	if($sQnameR1 ne $sQnameR2) {
		$iSingletons++;
		next;
	}
	if($sQnameR1 =~ /null/) {

	}
	if($sQnameR2 =~ /null/) {

	}
	ParseFlag($iFlagR1, \%hStatR1);
	ParseFlag($iFlagR2, \%hStatR2);
}

###############
# SUBROUTINES #
###############

sub ParseFlag {
	my ($iDec, $phStat) = @_;
	my ($iLbin, $iBin);
	$iLbin = unpack("B32", pack("N", $iDec));
	$iLbin =~ s/^0+(?=\d)//;    # otherwise you'll get leading zeros
	$iBin = sprintf("%012d", $iLbin);
	$iRbin = reverse($iBin);
        $phStat = {
			'paired' => substr($iRbin, 0, 1),
        	   	'propermap' => substr($iRbin, 1, 1),
        	   	'qunmapped' => substr($iRbin, 2, 1), 
        		'munmapped' => substr($iRbin, 3, 1), 
        		'qrev' => substr($iRbin, 4, 1), 
        		'mrev' => substr($iRbin, 5, 1), 
        		'firstpair' => substr($iRbin, 6, 1), 
        		'secondpair' => substr($iRbin, 7, 1), 
        		'scndryalign' => substr($iRbin, 8, 1), 
        		'failqual' => substr($iRbin, 9, 1), 
        		'pcrdup' => substr($iRbin, 10, 1), 
        		'supplealign' => substr($iRbin, 11, 1)
		  };
}

sub SortBam {
	my ($phCmdLineArgs, $sFile) = @_;
	my ($sCmd, $sOut);
	my $nExitCode;
	my $sOptions = "";	

	my $sSubName = (caller(0))[3];
	if(exists($phCmdLineArgs->{'samtools_params'})) {
		$sOptions = $phCmdLineArgs->{'samtools_params'};
	}
	$sOut = $phCmdLineArgs->{'output_dir'}."/".$sFile;
	$sCmd = $phCmdLineArgs->{'samtools_path'}." sort ".$sOptions." -O bam -n -o ".$sOut." -T /tmp/".$sFile." ".$phCmdLineArgs->{'input_file'};
	printLogMsg($DEBUG, "INFO : $sSubName :: Start sorting file $phCmdLineArgs->{'input_file'} to produce sorted BAM $sOut.\nINFO : $sSubName :: Command : $sCmd");
	$nExitCode = system($sCmd);
	if($nExitCode != 0) {
		printLogMsg($ERROR, "ERROR : $sSubName :: Sorting of BAM failed for $phCmdLineArgs->{'input_file'}/ with error. Check the stderr");
	} else {
		printLogMsg($DEBUG, "INFO : $sSubName :: Sorting of BAM succesfully completed in $phCmdLineArgs->{'output_dir'}/".$sFile);
	}
}


# Description   : Used to check the correctness of the command line arguments passed to the script. The script exits if required arguments are missing. 
# Parameters    : NA
# Returns       : NA
# Modifications :

sub checkCmdLineArgs {
	my ($phCmdLineArgs) = @_;
	my $sSubName = (caller(0))[3];
	if(exists($phCmdLineArgs->{'log'})) {
		open($fhLog, "> $phCmdLineArgs->{'log'}") or die "Could not open $phCmdLineArgs->{'log'} file for writing.Reason : $!\n"
	}
	my @aRequired = qw(input_file output_dir samtools_path);
        foreach my $sOption(@aRequired) {
                if(!defined($phCmdLineArgs->{$sOption})) {
                        printLogMsg($ERROR, "ERROR : $sSubName :: Required option $sOption not passed");
                }
        }
}

####################################################################################################################################################

# Description   : Used to handle logging of messages(errors and warnings) during the execution of the script
# Parameters    : level = can be ERROR, WARNING or INFO
#		  msg   = msg to be printed in the log file or to STDERR
# Returns       : NA
# Modifications : 

sub printLogMsg {
	my ($nLevel, $sMsg) = @_;
	if( $nLevel <= $DEBUG ) {
		print STDERR "$sMsg\n";
		print $fhLog "$sMsg\n" if(defined($fhLog));
		die "" if($nLevel == $ERROR);
	}	
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
