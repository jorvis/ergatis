#!/usr/bin/perl -w

eval 'exec /usr/bin/perl -w -S $0 ${1+"$@"}'
    if 0; # not running under some shell

#########################################################################################
#											#
# Name	      : download_sra.pl								#
# Version     : 1.0									#
# Project     :	LGT Seek Pipeline							#
# Description : Script to download files from SRA via FTP				#
# Author      : Sonia Agrawal								#
# Date        : January 23, 2014							#
#											#
#########################################################################################

use strict;
# Extended processing of command line options
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
# Prints usage message from embedded pod documentation
use Pod::Usage;
# Include packages here

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
my %hFTP = ();
my $phFTPConf = {};

my $fetch_metadata = 0;

################
# MAIN PROGRAM #
################
GetOptions(\%hCmdLineArgs,
	   'ftp|f=s',
	   'run_id|r=s',
	   'output_dir|o=s',
	   'username|u=s',
	   'password|p=s',
       'fetch_metadata|m=i',
	   'num_retry|n=i',
	   'log|l=s',
	   'help|h'
	  ) or pod2usage();

pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDERR} ) if ($hCmdLineArgs{'help'});

checkCmdLineArgs(\%hCmdLineArgs);

DownloadSRA(\%hCmdLineArgs);


###############
# SUBROUTINES #
###############

####################################################################################################################################################
# Description   : Used to download read files from NCBI SRA FTP server using wget
# Parameters    : phCmdLineArgs
#		  phCmdLineArgs - reference to hash of command line arguments passed to the perl script
# Returns       : NA
# Modifications :

sub DownloadSRA {
	my ($phCmdLineArgs) = @_;
	my ($sFile, $sCmd);
	my $nExitCode;
    my $sSubName = (caller(0))[3];
    foreach my $run_id ( split(/,\s*/, $phCmdLineArgs->{'run_id'}) ) {
        $sFile = CreateFilePath($run_id);

    # wget options :
    # ‘-r’ : Turn on recursive retrieving of directories.
    # ‘-nd’: Do not create a hierarchy of directories when retrieving recursively. With this option turned on, all files will get saved to the current directory, without clobbering (if a name shows up more than once, the filenames will get extensions ‘.n’).
    # ‘-c’ : Continue getting a partially-downloaded file.This is useful when you want to finish up a download started by a previous instance of Wget, incase a pipline is resumed.
    # ‘-N’ : Turn on time-stamping.
    # ‘-P prefix’ : Set directory prefix to prefix. The directory prefix is the directory where all other files and subdirectories will be saved to, i.e. the top of the retrieval tree.
    # ‘-t number’ :	Set number of retries to number. Specify 0 or ‘inf’ for infinite retrying. The default is to retry 20 times, with the exception of fatal errors like “connection refused” or “not found” (404), which are not retried.

        $sCmd = "wget -r -nd -c -N ";
        if(defined($phCmdLineArgs->{'num_retry'})) {
            $sCmd .= "-t $phCmdLineArgs->{'num_retry'} ";
        }
        $sCmd .= "-P $phCmdLineArgs->{'output_dir'} ftp://".$phCmdLineArgs->{'username'}.":".$phCmdLineArgs->{'password'}."@".$phCmdLineArgs->{'ftp'}."/".$sFile;

        printLogMsg($DEBUG, "INFO : $sSubName :: Start downloading $phCmdLineArgs->{'run_id'} in $phCmdLineArgs->{'output_dir'}.\nINFO : $sSubName :: Command : $sCmd");
        $nExitCode = system($sCmd);
    # wget returns 0 on success
        if($nExitCode == 0) {
            printLogMsg($DEBUG, "INFO : $sSubName :: $run_id downloaded succesfully in $phCmdLineArgs->{'output_dir'}");
        } else {
            printLogMsg($ERROR, "ERROR : $sSubName :: $run_id download failed with exit code $nExitCode. For details check STDERR file");
        }

        if ($fetch_metadata){
            $run_info_cmd = "wget -O $phCmdLineArgs->{'output_dir'} " . "/$phCmdLineArgs->{'run_id'}_info.csv http://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?save=efetch&db=sra&rettype=runinfo&term= $phCmdLineArgs->{'run_id'}";

            $nExitCode = system($run_info_cmd);
            # wget returns 0 on success
            if($nExitCode == 0) {
                printLogMsg($DEBUG, "INFO : $sSubName :: $run_id metadata downloaded succesfully in $phCmdLineArgs->{'output_dir'}");
            } else {
                printLogMsg($ERROR, "ERROR : $sSubName :: $run_id metadata download failed with exit code $nExitCode. For details check STDERR file");
            }
        }
	}
}

####################################################################################################################################################
# Description   : Used to create the file path on NCBI SRA FTP server using the SRX/SRR/SRS/SRP id
# Parameters    : sRunId
#		  sRunId - SRA compatible experiment id or run id or sample id or study id for the read files to be downloaded
# Returns       : sFile
#		  sFile - NCBI SRA FTP file path for the specified id
# Modifications :

sub CreateFilePath {
	my ($sRunId) = @_;
	my ($sFile, $sVol);
    	my $sSubName = (caller(0))[3];

	$sRunId =~ s/^\s+|\s+$//;
	$sFile = "/sra/sra-instant/reads/";
# Volume is the directory on FTP with first 6 characters from the SRA id
	$sVol = substr($sRunId, 0, 6);
	if($sRunId =~ /^SRX/) {
		$sFile .= "ByExp/sra/SRX/$sVol/$sRunId/*";
	} elsif($sRunId =~ /^SRR/) {
		$sFile .= "ByRun/sra/SRR/$sVol/$sRunId/*";
	} elsif($sRunId =~ /^SRS/) {
		$sFile .= "BySample/sra/SRS/$sVol/$sRunId/*";
	} elsif($sRunId =~ /^SRP/) {
		$sFile .= "ByStudy/sra/SRP/$sVol/$sRunId/*"
	} elsif($sRunId =~ /\/+/) {
# If user has specified the complete path of the file to be downloaded from FTP site
		$sFile = $sRunId;
	} else {
		printLogMsg($ERROR, "ERROR : $sSubName :: Unable to identify specified $sRunId id");
	}

	return($sFile);
}

####################################################################################################################################################
# Description   : Used to check the correctness of the command line arguments passed to the script. The script exits if required arguments are missing.
# Parameters    : phCmdLineArgs
#		  phCmdLineArgs - reference to hash of command line arguments passed to the perl script
# Returns       : NA
# Modifications :

sub checkCmdLineArgs {
	my ($phCmdLineArgs) = @_;
	my $sOption;
	my @aRequired = ();

	if(exists($phCmdLineArgs->{'log'})) {
		open($fhLog, "> $phCmdLineArgs->{'log'}") or die "Could not open $phCmdLineArgs->{'log'} file for writing.Reason : $!\n";
	}
	@aRequired = qw(ftp run_id output_dir);
        foreach $sOption(@aRequired) {
                if(!defined($phCmdLineArgs->{$sOption})) {
                        printLogMsg($ERROR,"ERROR! : Required option $sOption not passed");
                }
        }
	if(!exists($phCmdLineArgs->{'username'})) {
		$phCmdLineArgs->{'username'} = "anonymous";
	}
	if(!exists($phCmdLineArgs->{'password'})) {
		$phCmdLineArgs->{'password'} = "anonymous";
	}

    # Change to fetch_metadata param value if passed in.  Otherwise it's 0
    $fetch_metadata = $phCmdLineArgs->{'fetch_metadata'} if $phCmdLineArgs->{'fetch_metadata'};
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
		die "" if($nLevel == $ERROR);
	}
	print $fhLog "$sMsg\n" if(defined($fhLog));
}

__END__

#####################
# POD DOCUMENTATION #
#####################

=head1 NAME

download_sra.pl - Script to download read files from NCBI SRA FTP site

=head1 SYNOPSIS

# USAGE : perl download_sra.pl -r <SRA id> -f <FTP server> -o <path to output dir> [ -u <FTP server username> -p <FTP server password> -n <number of retries> -l <path to log file> ]

	parameters in [] are optional

=head1 OPTIONS

	-r <run_id>	:	NCBI SRA compatible 9-character id. Could be Study id (SRPXXXXXX), Experiment id (SRXXXXXXX), Run id (SRRXXXXXX) or Sample id (SRSXXXXXX). X stands for digit. Mandatory
	                Can separate multiple IDs with commas.

	-f <ftp>	:	Name of the NCBI FTP server. Currently it is ftp-trace.ncbi.nih.gov. Mandatory

	-o <output_dir>	:	Path to the output directory where the files will be downloaded. Mandatory

[
	-u <username>	:	Username for the FTP server. Default: anonymous. Optional

	-p <password>	:	Password for the FTP server. Default: anonymous. Optional

	-n <num_retry>	:	Number of retries to download a file. wget default is 20. Optional

    -m <fetch_metadata> :   If set to 1, will download CSV file of the SRA run info.  Optional

	-l <log>	: 	Path to log file. Optional

]

=head1 DESCRIPTION



=head1 INPUT

	NCBI SRA compatible valid 9-character ids, FTP server name and output directory where the files will be downloaded

=head1 OUTPUT

	Read files (.sra) from SRA

=head1 AUTHOR

	Sonia Agrawal
	Senior Bioinformatics Software Engineer
	Institute for Genome Sciences
	School of Medicine
	University of Maryland
	Baltimore, MD - 21201
	sagrawal@som.umaryland.edu

==cut
