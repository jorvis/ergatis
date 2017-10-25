#!/usr/bin/env perl -w

#####################################################################################
#
# Name	      : transform_mugsy_cog
# Version     :	1.0
# Project     :	CloVR Comparative Pipeline
# Description : Script to replace polypeptide ids in mugsycog file with GenBank gene 
#		ids
# Author      : Sonia Agrawal
# Date        : April 14, 2014
#
#####################################################################################

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
my %hMap = ();
# Log file handle;
my $fhLog;
my ($ERROR, $WARN, $DEBUG) = (1, 2, 3);

################
# MAIN PROGRAM #
################
GetOptions(\%hCmdLineArgs,
	   'mugsycog|c=s',
	   'map_file|m=s',
	   'output_file|o=s',
	   'log|l=s', 
	   'help|h'
	  ) or pod2usage();

pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDERR} ) if ($hCmdLineArgs{'help'});

checkCmdLineArgs();

readMapFile($hCmdLineArgs{'map_file'}, \%hMap);

replacePolyIds($hCmdLineArgs{'mugsycog'}, $hCmdLineArgs{'output_file'}, \%hMap);


###############
# SUBROUTINES #
###############

# Description   : Used to check the correctness of the command line arguments passed to the script. The script exits if required arguments are missing. 
# Parameters    : NA
# Returns       : NA
# Modifications :

sub readMapFile {
	my ($sFile, $phMap) = @_;
	my ($sLine, $sMol, $sGene);
	my @aCols = (); 

	open(FR, "< $sFile") or printLogMsg($ERROR, "ERROR! : Could not open $sFile file for reading. Reason : $!");
	
	while($sLine=<FR>) {
		chomp($sLine);
		next if($sLine =~ /^\s*$/);
		@aCols = split(/\t/,$sLine);
		($sMol, $sGene) = split(/\|\|\|/, $aCols[0], 2);
		$phMap->{$aCols[5]} = $sMol.".".$sGene;
	}

	close(FR);
}

####################################################################################################################################################
# Description   : Used to check the correctness of the command line arguments passed to the script. The script exits if required arguments are missing. 
# Parameters    : NA
# Returns       : NA
# Modifications :

sub replacePolyIds {
	my ($sCogFile, $sOutFile, $phMap) = @_;
	my ($sLine);
	
	open(FR, "< $sCogFile") or printLogMsg($ERROR, "ERROR! : Could not open $sCogFile file for reading. Reason : $!");
	open(FW, "> $sOutFile") or printLogMsg($ERROR, "ERROR! : Could not open $sOutFile file for writing. Reason : $!");
	
	while($sLine=<FR>) {
		chomp($sLine);
		if($sLine =~ /^COG/) {
			print FW "$sLine\n";
		} elsif($sLine =~ /^(\s+)(\S+)/) {
			if(exists($phMap->{$2})) {
				print FW $1.$phMap->{$2}."\n";
			} else {
				print FW "$sLine\n";
				printLogMsg($WARN, "WARNING : Unable to find mappping for $2. Please check the input $sCogFile and mapping files");
			}
		} else {
			print FW "$sLine\n";
		}
	}			
	close(FR);
	close(FW);
}

####################################################################################################################################################
# Description   : Used to check the correctness of the command line arguments passed to the script. The script exits if required arguments are missing. 
# Parameters    : NA
# Returns       : NA
# Modifications :

sub checkCmdLineArgs {
	if(exists($hCmdLineArgs{'log'})) {
		open($fhLog, "> $hCmdLineArgs{'log'}") or die "Could not open $hCmdLineArgs{'log'} file for writing.Reason : $!\n"
	}
	my @aRequired = qw(mugsycog map_file output_file);
        foreach my $sOption(@aRequired) {
                if(!defined($hCmdLineArgs{$sOption})) {
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
