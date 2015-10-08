#!/usr/local/bin/perl -w

#########################################################################################
#											#
# Name	      : sra2fastq.pl								#
# Version     : 1.0									#
# Project     :	LGT Pipeline								#
# Description : Script to convert downloaded SRA files from SRA into FASTQ format	#
# Author      : Sonia Agrawal								#
# Date        : August 31, 2015								#
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
my $logfh;
my ($ERROR, $WARN, $DEBUG) = (1, 2, 3);

################
# MAIN PROGRAM #
################
GetOptions(\%hCmdLineArgs,
	   'input_file|i=s',
	   'sratoolkit|s=s',
	   'output_dir|o=s',
	   'paired|p=i',
	   'log|l=s', 
	   'help|h'
	  ) or pod2usage();

pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDERR} ) if ($hCmdLineArgs{'help'});

checkCmdLineArgs(\%hCmdLineArgs);

ConvertSRS(\%hCmdLineArgs);

###############
# SUBROUTINES #
###############

####################################################################################################################################################
# Description   : Used to convert .sra files downloaded from SRA FTP server to FASTQ format
# Parameters    : phCmdLineArgs 
#		  phCmdLineArgs - reference to hash of command line arguments passed to the perl script
# Returns       : NA
# Modifications :

sub ConvertSRS {
	my ($phCmdLineArgs) = @_;
	my ($sCmd);
	my $nExitCode;
    	my $sSubName = (caller(0))[3];
# Convert to fastq format
	$sCmd = $phCmdLineArgs->{'sratoolkit'}."/bin/fastq-dump --split-3 -O ".$phCmdLineArgs->{'output_dir'}." ".$phCmdLineArgs->{'input_file'};
	printLogMsg($DEBUG, "INFO : $sSubName :: Start converting $phCmdLineArgs->{'input_file'} to FASTQ files in $phCmdLineArgs->{'output_dir'}.\nINFO : $sSubName :: Command : $sCmd");
	$nExitCode = system($sCmd);
	#fastq-dump returns 0 on success
	if($nExitCode != 0) {
		printLogMsg($ERROR, "ERROR : $sSubName :: $phCmdLineArgs->{'input_file'} conversion failed with error");
	} else {
		printLogMsg($DEBUG, "INFO : $sSubName :: $phCmdLineArgs->{'input_file'} SRA to FASTQ conversion succesfully completed in $phCmdLineArgs->{'output_dir'}");
	}

# Convert to fasta format
    $sCmd = $phCmdLineArgs->{'sratoolkit'}."/bin/fastq-dump --fasta 60 -O ".$phCmdLineArgs->{'output_dir'}." ".$phCmdLineArgs->{'input_file'} .;
    printLogMsg($DEBUG, "INFO : $sSubName :: Start converting $phCmdLineArgs->{'input_file'} to FASTA files in $phCmdLineArgs->{'output_dir'}.\nINFO : $sSubName :: Command : $sCmd");
    $nExitCode = system($sCmd);
    #fastq-dump returns 0 on success
    if($nExitCode != 0) {
        printLogMsg($ERROR, "ERROR : $sSubName :: $phCmdLineArgs->{'input_file'} conversion failed with error");
    } else {
        printLogMsg($DEBUG, "INFO : $sSubName :: $phCmdLineArgs->{'input_file'} SRA to FASTA conversion succesfully completed in $phCmdLineArgs->{'output_dir'}");

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
		open($logfh, "> $phCmdLineArgs->{'log'}") or die "Could not open $phCmdLineArgs->{'log'} file for writing.Reason : $!\n";
	}
	@aRequired = qw(input_file sratoolkit output_dir);
        foreach $sOption(@aRequired) {
                if(!defined($phCmdLineArgs->{$sOption})) {
                        printLogMsg($ERROR,"ERROR! : Required option $sOption not passed");
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
	my ($level, $msg) = @_;
	if( $level <= $DEBUG ) {
		print STDERR "$msg\n";
		print $logfh "$msg\n" if(defined($logfh));
		die "" if($level == $ERROR);
	}	
}

__END__

#####################
# POD DOCUMENTATION #
#####################

=head1 NAME

sra2fastq.pl - Script to convert .sra files downloaded from NCBI SRA FTP site to FASTQ format

=head1 SYNOPSIS

# USAGE : perl sra2fastq.pl -i <path to input file> -s <path to SRA took kit fastq-dump binary> -o <path to output dir> [ -l <path to log file> ]

	parameters in [] are optional

=head1 OPTIONS

	-i <input_file>	:	Path to the input SRA file to be converted to FASTQ. Mandatory
	
	-s <sratoolkit>	:	Path to sratoolkit binary for fastq-dump utility. Mandatory

	-o <output_dir>	:	Path to the output directory where the files will be downloaded. Mandatory

[	

	-l <log>	: 	Path to log file. Optional
] 

=head1 DESCRIPTION



=head1 INPUT

	Path to SRA file downloaded from NCBI SRA, path to sratoolkit fastq-dump binary and output directory where the files will be downloaded

=head1 OUTPUT

	Read files (.fastq) 
	1. First biological reads are placed in files *_1.fastq and 
        2. *_2.fastq.
	3. If only one biological read is present it is placed in *.fastq

=head1 AUTHOR

	Sonia Agrawal
	Senior Bioinformatics Software Engineer
	Institute for Genome Sciences
	School of Medicine
	University of Maryland
	Baltimore, MD - 21201
	sagrawal@som.umaryland.edu

==cut
