#!/usr/local/bin/perl -w

#########################################################################################
#											#
# Name	      : tbl2asn.pl								#
# Version     : 1.0									#
# Project     :	GenBank Submisison Pipeline						#
# Description : Script to run NCBI tbl2asn 						#
# Author      : Sonia Agrawal								#
# Date        : February 10, 2014							#
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
my ($sCmd, $sDiscrep, $sSource, $sSqnFile, $sTbl, $sFasta, $sSbt);
my @aMeta = ();

################
# MAIN PROGRAM #
################
GetOptions(\%hCmdLineArgs,
	   'input_file|i=s',
	   'input_dir|d=s',
	   'output_dir|o=s',
	   'utility_path|p=s',
	   'opts|t=s',
	   'log|l=s', 
	   'help|h'
	  ) or pod2usage();

pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDERR} ) if ($hCmdLineArgs{'help'});

checkCmdLineArgs(\%hCmdLineArgs);

readInput($hCmdLineArgs{'input_file'}, \@aMeta);

$sDiscrep = $hCmdLineArgs{'output_dir'}."/".$aMeta[0]."_discrep.txt";

if(-e $hCmdLineArgs{'input_dir'} && -d $hCmdLineArgs{'input_dir'}) {
	$sTbl = $hCmdLineArgs{'input_dir'}."/".$aMeta[0].".tbl";
	$sFasta = $hCmdLineArgs{'input_dir'}."/".$aMeta[0].".fsa";
	$sSbt = $hCmdLineArgs{'input_dir'}."/".$aMeta[0].".sbt";
	if(-e $sFasta && -e $sSbt) {
		$sSource = "[gcode=11][host=$aMeta[8]][country=$aMeta[10]][collection-date=$aMeta[9]][organism=$aMeta[5]][strain=$aMeta[6]][serotype=$aMeta[7]][isolation-source=$aMeta[19]][tech=wgs]";
		if((-e $hCmdLineArgs{'utility_path'}) && (-x $hCmdLineArgs{'utility_path'})) {
			$sCmd = "$hCmdLineArgs{'utility_path'} -p $hCmdLineArgs{'input_dir'} -t $sSbt -r $hCmdLineArgs{'output_dir'} -a s -V vb -X C -Z $sDiscrep -j \"$sSource\"";
			if(defined($hCmdLineArgs{'opts'})) {
				$sCmd = $sCmd." ".$hCmdLineArgs{'opts'};
			}
		} else {
			printLogMsg($ERROR, "ERROR : Path to tbl2asn script $hCmdLineArgs{'utility_path'} is invalid or the script is not executable.");
		}	
	} else {
		printLogMsg($ERROR, "ERROR : Specified input directory $hCmdLineArgs{'input_dir'} is missing .tbl, .fsa and/or .sbt files required for running tbl2asn.");
	}	
} else {
	printLogMsg($ERROR, "ERROR : Specified input directory $hCmdLineArgs{'input_dir'} does not exist or is not a valid directory.");
}



printLogMsg($DEBUG, "INFO : Changing directory to output directory $hCmdLineArgs{'output_dir'}");
my $sRet = chdir($hCmdLineArgs{'output_dir'});

printLogMsg($DEBUG, "INFO : Executing tbl2asn command :: $sCmd");
system($sCmd);

$sSqnFile = $hCmdLineArgs{'output_dir'}."/".$aMeta[0].".sqn";
if((-e $sSqnFile) && (-s $sSqnFile)) {
	printLogMsg($DEBUG, "INFO : tbl2asn ran successfully created files in output directory $hCmdLineArgs{'output_dir'}");
} else {
	printLogMsg($ERROR, "ERROR : Error executing tbl2asn. Required files not created!")
}

###############
# SUBROUTINES #
###############

####################################################################################################################################################
# Description   : Used to read input meta data file
# Parameters    : sFile, paMeta
#		  sFile - Path to input meta data file containing the information as described below
#		  paMeta - pointer to array to hold meta data
# Returns       : NA
# Assumptions	: The input meta data file sFile contains only a single row with information about a single database. This script is designed to be
#		  run interatively within an ergatis component. So it processes one DB at a time. The ergatis component itself could process 
#		  several databases.
# Modifications :

sub createDbDir {
        my ($sOutDir, $sDb) = @_; 
        my $sSubName = (caller(0))[3];
        my $sDbDir;

        $sDbDir = $sOutDir."/".$sDb;
        if((-e $sDbDir) && (-d $sDbDir)) {
                return(1);
        } else {
                if((-e $sOutDir) && (-d $sOutDir)) {
                        if(mkdir $sDbDir) {
                                return(1);
                        } else {
                                printLogMsg($ERROR, "ERROR : $sSubName :: Unable to create database specific directory within output directory $sOutDir. Reason : $!");
                        }   
                            
                } else {
                        printLogMsg($ERROR, "ERROR : $sSubName :: Provided output directory $sOutDir does not exist.");
                }   
        }
}

####################################################################################################################################################
# Description   : Used to read input meta data file
# Parameters    : sFile, paMeta
#		  sFile - Path to input meta data file containing the information as described below
#		  paMeta - pointer to array to hold meta data
# Returns       : NA
# Assumptions	: The input meta data file sFile contains only a single row with information about a single database. This script is designed to be
#		  run interatively within an ergatis component. So it processes one DB at a time. The ergatis component itself could process 
#		  several databases.
# Modifications :

sub readInput {
	my ($sFile, $paMeta) = @_;
	my $nI;
	my $sLine;
	my $fhRead;
    	my $sSubName = (caller(0))[3];
	my @aMandatory = (0, 5, 6);

	open($fhRead, "< $sFile") or printLogMsg($ERROR, "ERROR : $sSubName :: Could not open $sFile file for reading. Reason : $!");
	while($sLine=<$fhRead>) {
		chomp($sLine);
		next if($sLine =~ /^#/);
		next if($sLine =~ /^\s+$/);
		@{$paMeta} = split(/\t/, $sLine);
		# @meta : This script needs columns 0 and 5-10 
#		[0] = Db name
#		[1] = NCBI locus tag
#		[2] = Path to rules file
#		[3] = Path to gene symbols file
#		[4] = Bioproject Id
#		[5] = Organism name
#		[6] = strain name
#		[7] = Organism serotype
#		[8] = Host
#		[9] = Date
#		[10]= Country
#		[11]= Assembly method
#		[12]= Genome Coverage
#		[13]= Sequencing platform
#		[14]= Contact person's name. Last name\sFirst name
#		[15]= Contact person's email address
#		[16]= Comma-separated author list. Each name should be Last name\sFirst name\sMiddle Initial
#		[17]= Title of the publication
#		[18]= Path to EC number rules file for EC numbers to be deleted
#		[19]= Source of isolation of the sample
		foreach $nI (@aMandatory) {
			if(length($paMeta->[$nI]) == 0 ) {
				printLogMsg($ERROR, "ERROR : $sSubName :: $nI column value missing in $sLine.");
			} 
		}
		if($paMeta->[19] eq "NA") {
			$paMeta->[19] = "";
		}
	}
	close($fhRead);	
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
    	my $sSubName = (caller(0))[3];

	if(exists($phCmdLineArgs->{'log'})) {
		open($logfh, "> $phCmdLineArgs->{'log'}") or die "ERROR : $sSubName :: Could not open $phCmdLineArgs->{'log'} file for writing.Reason : $!\n";
	}
	@aRequired = qw(input_file input_dir output_dir utility_path);
        foreach $sOption(@aRequired) {
                if(!defined($phCmdLineArgs->{$sOption})) {
                        printLogMsg($ERROR,"ERROR : $sSubName :: Required option $sOption not passed");
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

fix_gene_symbols.pl - Script to remove upper case, duplicate gene symbols as well replace existing gene symbols in DB with new ones provided in the rules file

=head1 SYNOPSIS

# USAGE : perl fix_gene_symbols.pl -i <meta data file> [ -o <path to output dir>  -u <database username> -p <database password> -s <database server> -t <path to .tbl file> -l <path to log file> ]

	parameters in [] are optional

=head1 OPTIONS

	-i <input_file>	:	Path to input meta data file containing one row for a database. Mandatory
[	
	-o <output_dir>	:	Path to the output directory where the corrected tbl file will be created. Mandatory if change required in only tbl file

	-u <username>	:	Username for the database server. Mandatory if change required in the database only.

	-p <password>	:	Password for the database server. Mandatory if change required in the database only.

	-s <host>	:	Name of the database server. Mandatory if change required in the database only.

	-t <tbl_file>	:	Path to the original .tbl file to be corrected for gene symbols. Mandatory if change required in only tbl file

	-l <log>	: 	Path to log file. Optional

] 

=head1 DESCRIPTION



=head1 INPUT

	Required input is the meta data file with the following columns in a row for a database:
		[0] = Db name
		[1] = NCBI locus tag
		[2] = Path to rules file
		[3] = Path to gene symbols file
		[4] = Bioproject Id
		[5] = Organism name
		[6] = strain name
		[7] = Organism serotype
		[8] = Host
		[9] = Date
		[10]= Country
		[11]= Assembly method
		[12]= Genome Coverage
		[13]= Sequencing platform
		[14]= Contact person's name. Last name\sFirst name
		[15]= Contact person's email address
		[16]= Comma-separated author list. Each name should be Last name\sFirst name\sMiddle Initial
		[17]= Title of the publication
	For this script only column [0] and [3] needs to be have valid values.

	If .tbl file needs to be changed, then provide the original tbl file dumped from manatee database. If directly database needs to changed 
	then pass the required database parameters to the script. Either .tbl file or database parameters needs to be present but not both.

=head1 OUTPUT

	Corrected .tbl file if original tbl file was passed. Otherwise, gene symbol changes happen in the database provided in the input meta data file.

=head1 AUTHOR

	Sonia Agrawal
	Senior Bioinformatics Software Engineer
	Institute for Genome Sciences
	School of Medicine
	University of Maryland
	Baltimore, MD - 21201
	sagrawal@som.umaryland.edu

==cut
