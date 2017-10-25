#!/usr/bin/env perl -w

#########################################################################################
#											#
# Name	      : fix_gbk									#
# Version     : 1.0									#
# Project     : Clovr Microbe pipeline							#
# Description : Script to add accession ids to gbk and fix taxonomy			#
# Author      : Sonia Agrawal								#
# Date        : June 30, 2014								#		
#											#
#########################################################################################

use strict;
# Extended processing of command line options
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
# Prints usage message from embedded pod documentation
use Pod::Usage;
# Include packages here 
use Ergatis::IdGenerator;

#############
# CONSTANTS #
#############


###########
# GLOBALS #
###########
my %cmdLineArgs = ();
# Log file handle;
my $logfh;
my ($ERROR, $WARN, $DEBUG) = (1, 2, 3);
my $oIdCreator;
my ($sLine, $sId, $sLocus);
my $iFlag = 0;

################
# MAIN PROGRAM #
################
GetOptions(\%cmdLineArgs,
	   'input_file|i=s',
	   'id_repository|r=s',
	   'output_file|o=s',
	   'log|l=s', 
	   'help|h'
	  ) or pod2usage();

pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDERR} ) if ($cmdLineArgs{'help'});

checkCmdLineArgs();

$oIdCreator = new Ergatis::IdGenerator('id_repository' => $cmdLineArgs{'id_repository'});

open(FR, "< $cmdLineArgs{'input_file'}") or printLogMsg($ERROR, "ERROR!: Could not open $cmdLineArgs{'input_file'} file for reading. Reason : $!");

open(FW, "> $cmdLineArgs{'output_file'}") or printLogMsg($ERROR, "ERROR!: Could not open $cmdLineArgs{'output_file'} file for writing. Reason : $!");


while($sLine = <FR>) {
	chomp($sLine);
	if($sLine =~ /^LOCUS/) {
		$sLocus = $sLine;
		next;
	}
	if(length($sLocus) > 0) {
		if($sLine !~ /^DEFINITION/) {
			$sLine = $sLocus." ".$sLine;
		} else {
			print FW $sLocus."\n";
		}
		$sLocus = "";			
	}
	if($sLine =~ /^ACCESSION\s*/) {
		$sId = $oIdCreator->next_id('type' => 'accession','project' => 'test');
		$sId =~ s/test\.accession\.//;
		print FW "ACCESSION". ' ' x 3 . "$sId\n";
		next;
	}
	if($sLine =~ /^\s+ORGANISM/) {
		$iFlag = 1;
		print FW $sLine."\n";
		next;
	}
	if($iFlag) {
		print FW ' ' x 12 . "Bacteria; Unclassified.\n";
		$iFlag = 0;
		next;
	}
	print FW $sLine."\n";
}

close(FR);
close(FW);

###############
# SUBROUTINES #
###############

# Description   : Used to check the correctness of the command line arguments passed to the script. The script exits if required arguments are missing. 
# Parameters    : NA
# Returns       : NA
# Modifications :

sub checkCmdLineArgs {
	if(exists($cmdLineArgs{'log'})) {
		open($logfh, "> $cmdLineArgs{'log'}") or die "Could not open $cmdLineArgs{'log'} file for writing.Reason : $!\n"
	}
	my @required = qw(input_file id_repository output_file);
        foreach my $option(@required) {
                if(!defined($cmdLineArgs{$option})) {
                        printLogMsg($ERROR,"ERROR! : Required option $option not passed");
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

fix_gbk.pl - Script to add accession numbers and fix taxonomy information in newly created GenBank file thorough the clovr microbe pipeline

=head1 SYNOPSIS

# USAGE : perl fix_gbk.pl --i <Path to input GenBank file> --o <Path to output GenBank file> --r <Path to Ergatis id repository> [--l <Path to log file>]

	parameters in [] are optional

=head1 OPTIONS



=head1 DESCRIPTION



=head1 INPUT



=head1 OUTPUT



=head1 AUTHOR

	Sonia Agrawal
	Senior Bioinformatics Software Engineer
	Institute for Genome Sciences
	School of Medicine
	University of Maryland
	Baltimore, MD - 21201
	sagrawal@som.umaryland.edu

==cut
