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
use File::Basename;

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

my $data_dir_present = 0;
my $output_dir;

################
# MAIN PROGRAM #
################
GetOptions(\%hCmdLineArgs,
	   'input_file|i=s',
	   'sratoolkit|s=s',
	   'output_dir|o=s',
	   'data_dir|d=s',
	   'log|l=s',
	   'help|h'
	  ) or pod2usage();

pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDERR} ) if ($hCmdLineArgs{'help'});

checkCmdLineArgs(\%hCmdLineArgs);
$data_dir_present = is_data_dir_present(\%hCmdLineArgs);
ConvertSRS(\%hCmdLineArgs);

###############
# SUBROUTINES #
###############

sub is_data_dir_present{
	my $opts = shift;
	$output_dir = $opts->{'output_dir'};
	if ($opts->{'data_dir'}){
		$output_dir = $opts->{'data_dir'};
		return 1;
	}
	return 0;
}

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
	$sCmd = $phCmdLineArgs->{'sratoolkit'}."/bin/fastq-dump --split-3 -O " . $output_dir . " ".$phCmdLineArgs->{'input_file'};
	printLogMsg($DEBUG, "INFO : $sSubName :: Start converting $phCmdLineArgs->{'input_file'} to FASTQ files in $output_dir.\nINFO : $sSubName :: Command : $sCmd");
	$nExitCode = system($sCmd);
	#fastq-dump returns 0 on success
	if($nExitCode != 0) {
		printLogMsg($ERROR, "ERROR : $sSubName :: $phCmdLineArgs->{'input_file'} conversion failed with error");
	} else {
		printLogMsg($DEBUG, "INFO : $sSubName :: $phCmdLineArgs->{'input_file'} SRA to FASTQ conversion succesfully completed in $output_dir");
	}

	# Creating a blank file using the SRA ID.
	# Reason:  lgt_bwa component accepts an input_directory of fastq files, so iterating over
	# 			a fastq list with paired-end fastq files results in two groups iterating over
	# 			the same directory.  This ensures the directory is iterated over just once in
	# 			that component.
	my ($base, $dir, $ext) = fileparse($phCmdLineArgs->{'input_file'}, qr/\.[^.]*/);

	# If we need to symlink our contents over from the 'data directory', find the newly created fastq files and perform the symlink.  I'm having trouble figuring how how to best do this in an --exec so I'm writing in a script and executing that.
	if ($data_dir_present) {
		open SCRIPT, $phCmdLineArgs->{'output_dir'} . "ln_cmd.sh" || die "Cannot create script for symlinking: $!\n";
		(my $script = <<"END_SCRIPT") =~ s/^ {8}//gm;

		#!/usr/bin/bash

		base=basename $1
		ln -s $1 $phCmdLineArgs->{'output_dir'}\/${base}

END_SCRIPT
		print SCRIPT $script;
		close SCRIPT;
		my $symCmd = "find $output_dir --name ${base}*\.fastq --exec " . $phCmdLineArgs->{'output_dir'} . "/ln_cmd.sh {} \;";
		$nExitCode = system($symCmd);
		#fastq-dump returns 0 on success
		if($nExitCode != 0) {
			printLogMsg($ERROR, "ERROR : $sSubName :: Unable to symlink output FASTQ files");
		} else {
			printLogMsg($DEBUG, "INFO : $sSubName :: Symlinking FASTQ files to " . $phCmdLineArgs->{'output_dir'} . " successful");
		}
	}

	my $blank_file = $phCmdLineArgs->{'output_dir'} . "/" . $base . ".blank";
	`touch $blank_file`;
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
		die "" if($level == $ERROR);
	}
	print $logfh "$msg\n" if(defined($logfh));
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
	-d <data_dir>	:	Path to a directory to dump the converted fastq files.  The data will be symlinked to the <output_dir> to make it easier to manage in a pipeline.  Optional

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
