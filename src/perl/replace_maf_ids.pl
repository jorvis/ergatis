#!/usr/local/bin/perl -w

#################################################################################
#										#
# Name	      :	replace_maf_ids.pl						#
# Version     : 1.0								#
# Project     : CloVR comparative pipeline					#
# Description : Replace checksum ids with organism names in MAF file		# 	
# Author      : Sonia Agrawal							#
# Date        : April 24, 2013							#
#										#
#################################################################################

use strict;
# Extended processing of command line options
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
# Prints usage message from embedded pod documentation
use Pod::Usage;
# Include packages here 


###########
# GLOBALS #
###########
my %cmdLineArgs = ();
# Log file handle;
my $logfh;
my $id_map = {};
my $line;
my ($ERROR, $WARN, $DEBUG) = (1, 2, 3);

################
# MAIN PROGRAM #
################
GetOptions(\%cmdLineArgs,
	   'id_file|i=s',
	   'maf|m=s',
	   'output_file|o=s',
	   'log|l=s', 
	   'help|h'
	  ) or pod2usage();

pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDERR} ) if ($cmdLineArgs{'help'});

checkCmdLineArgs();

$id_map = readIdFile($cmdLineArgs{'id_file'});

open(FR, "< $cmdLineArgs{'maf'}") or printLogMsg($ERROR, "ERROR! : Could not open file $cmdLineArgs{'maf'} for reading. Reason: $!");
open(FW, "> $cmdLineArgs{'output_file'}") or printLogMsg($ERROR, "ERROR! : Could not open file $cmdLineArgs{'output_file'} for writing. Reason: $!");

while($line = <FR>) {
	chomp($line);
	if($line =~ /^s\s+(\w+)\.(\S+)(.+)/) {
		if(exists($id_map->{$1})) {
			print FW "s $id_map->{$1}.$2"."$3\n";
		} else {
			printLogMsg($ERROR, "ERROR! : Could not map organism name for $1.$2");
		}	
	} else {
		print FW "$line\n";
	}
}
close(FR);
close(FW);

###############
# SUBROUTINES #
###############

# Description   : Used to read in id file and create a hash. 
# Parameters    : Path to id file
# Returns       : Reference to hash of ids and organism names.Key=checksum, Value=organism name
# Modifications :

sub readIdFile {
	my ($filename) = @_;
	my ($line, $checksum, $id);
	my %mapping = ();
	open(FR, "< $filename") or printLogMsg($ERROR, "ERROR! : Could not open file $filename for reading. Reason: $!");
	while($line = <FR>) {
		chomp($line);
		next if($line =~ /^\s+$/);
		($id, $checksum) = split(/\t/, $line, 2);
		$mapping{$checksum} = $id;
	}
	close(FR);
	return(\%mapping);
}

####################################################################################################################################################

# Description   : Used to check the correctness of the command line arguments passed to the script. The script exits if required arguments are missing. 
# Parameters    : NA
# Returns       : NA
# Modifications :

sub checkCmdLineArgs {
	if(exists($cmdLineArgs{'log'})) {
		open($logfh, "> $cmdLineArgs{'log'}") or die "Could not open $cmdLineArgs{'log'} file for writing.Reason : $!\n"
	}
	my @required = qw();
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

replace_maf_ids.pl - Script to replace checksum ids from MAF file with organism names using a mapping file produced earlier in prepmugsy step.

=head1 SYNOPSIS

USAGE : perl replace_maf_ids.pl -i <id_file> -m <maf> -o <output_file> [-l <log_file>]

	parameters in [] are optional

=head1 OPTIONS
	
	--id_file	=	Path to a two-column tab delimited map file containing organims name\tchecksum

	--maf		=	Path to Multiple Alignment File (MAF) produced by Mugsy

	--output_file	=	Path to output MAF with replaced ids

	--log		=	Path to log file for logging errors, warnings 	


=head1 DESCRIPTION

	During the mugsyprep step in comparative genomics pipeline, organism names are replaced with checksums to avoid errors due to special characters in organism names. After runnig mugsyalign these 
	ids have to replaced back so that downstream steps will show logical organism names and not checksums.

=head1 INPUT

	1. map_id.txt - Two-column tab delimited file where first column is organism name from GenBank file and second column is checksum generated during mugsyprep step. Eg
		Neisseria_meningitidis_MC58     28975cbc1be5a66a45f3bbf3b75a7a48
		Neisseria_meningitidis_Z2491    dafa738d29791390e27a31790cfe54b1
 
	2. MAF - Multiple Alignment File produced by MUGSY containing whole genome multiple alignments.

=head1 OUTPUT

	1. output_file.maf - Path to output MAF file which will be same as input MAF with checksum ids replaced with organism name.

=head1 AUTHOR

	Sonia Agrawal
	Bioinformatics Software Engineer II
	Institute for Genome Sciences
	School of Medicine
	University of Maryland
	Baltimore, MD - 21201
	sagrawal@som.umaryland.edu

==cut
