#!/usr/local/bin/perl -w

#########################################################################################
#											#
# Name	      : prepare_sbt_cmt.pl							#
# Version     : 1.0									#
# Project     :	GenBank Submisison Pipeline						#
# Description : Script to create .sbt and .cmt files required for GenBank submission	#
# Author      : Sonia Agrawal								#
# Date        : January 31, 2014							#
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

# This is the URI link to pass params to in order to generate a .sbt file
my $SUBMISSION_TEMPLATE = 'https://submit.ncbi.nlm.nih.gov/genbank/template/submission/';
#my $SUBMISSION_TEMPLATE = 'http://www.ncbi.nlm.nih.gov/WebSub/template.cgi';

###########
# GLOBALS #
###########
my %hCmdLineArgs = ();
# Log file handle;
my $logfh;
my ($ERROR, $WARN, $DEBUG) = (1, 2, 3);
my %hFTP = ();
my @aMeta = ();

################
# MAIN PROGRAM #
################
GetOptions(\%hCmdLineArgs,
	   'input_file|i=s',
	   'output_dir|o=s',
	   'log|l=s', 
	   'help|h'
	  ) or pod2usage();

pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDERR} ) if ($hCmdLineArgs{'help'});

checkCmdLineArgs(\%hCmdLineArgs);

readInput($hCmdLineArgs{'input_file'}, \@aMeta);

prepareSbt(\@aMeta, $hCmdLineArgs{'output_dir'});

prepareCmt(\@aMeta, $hCmdLineArgs{'output_dir'});

###############
# SUBROUTINES #
###############

####################################################################################################################################################
# Description   : Used to create the file path on NCBI SRA FTP server using the SRX/SRR/SRS/SRP id  
# Parameters    : sRunId
#		  sRunId - SRA compatible experiment id or run id or sample id or study id for the read files to be downloaded
# Returns       : sFile
#		  sFile - NCBI SRA FTP file path for the specified id
# Modifications :

sub prepareSbt {
	my ($paMeta, $sOutDir) = @_;
	my $sSubName = (caller(0))[3];
	my $nI;
	my ($sCmd, $sSbtFile);
	my $fhRead;

	$sSbtFile = $sOutDir."/".$paMeta->[0].".sbt";
# Authors section
	my @aAuthors = ();
	my ($sLast, $sFirst, $sMiddle, $sAuthorForm) = "";
	my $nCnt;
	@aAuthors = split(/\,/, $paMeta->[16]);	
	for($nI=0; $nI < @aAuthors; $nI++) {
		$sLast = "";
		$sFirst = "";
		$aAuthors[$nI] =~ s/^\s+|\s+$//;
		($sLast, $sFirst, $sMiddle) = split(/\s+/, $aAuthors[$nI], 3);
		if(!defined($sMiddle)) {
			$sMiddle = "";
		} 	
		$sLast =~ s/^\s+|\s+$//;
		$sFirst =~ s/^\s+|\s+$//;
		$sMiddle =~ s/^\s+|\s+$// if(length($sMiddle) > 0);
		$nCnt = $nI + 1;
		$sAuthorForm .= " -F 'sequence_author-$nCnt-first_name=$sFirst' -F 'sequence_author-$nCnt-middle_initial=$sMiddle' -F 'sequence_author-$nCnt-last_name=$sLast' -F 'sequence_author-$nCnt-suffix='";
	}

# Contact information section
	my ($sCLast, $sCFirst);
	($sCLast, $sCFirst) = split(/\s+/, $paMeta->[14], 2);
	
	if (!$sCLast || !$sCFirst || !$paMeta->[15]) {
		printLogMsg($ERROR, "Author name (last, first) and author e-mail are required in metadata file: $!");
	}

	# Initializing hash of contact infomation with default values if user-provided information is not present
	my %contact_info;
	if (! $paMeta->[20] || ! -e $paMeta->[20]) {
		%contact_info = (
		'city' => 'Baltimore',
		'country' => 'USA',
		'department' => 'Institute for Genome Sciences',
		'fax' => 'none',
		'institution' => 'University of Maryland School of Medicine',
		'phone' => '410-706-1401',
		'state' => 'MD',
		'street' => 'BioPark II, 801 W. Baltimore St., Suite 619',
		'zip' => '21201'
		);
	} else {
		# Open the file of contacts, parse the entry for the author we need and write to hash
		open(CONTACT, $paMeta->[20]) || printLogMsg($ERROR, "Could not open file $paMeta->[20] for reading: $!");
		while (<CONTACT>) {
			chomp;
			next unless (/^$paMeta->[14]/);
			my @fields = split(/\t/);
			for (my $f=0; $f <=$#fields; $f++){
				# give any undefined fields a 'none' value
				$fields[$f] = 'none' if (! $fields[$f]);
			}
			# Populate the contact hash by iterating through the fields array
			my $field_index = 1;
			foreach my $field qw(institution department city state country street fax phone zip) {
				$contact_info{$field} = $fields[$field_index++];
			}
			foreach my $required qw(institution street city country phone) {
				printLogMsg($ERROR, "Field [$required] is required to have information in the contacts file: $!") if $contact_info{$required} eq 'none'; 
			}
			last;
		}
		close CONTACT;
	}

# Running NCBI template.cgi using cUrl
	$sCmd = "curl -F 'first_name=$sCFirst' -F 'last_name=$sCLast' -F 'department=$contact_info{department}' -F 'institution=$contact_info{institution}' -F 'street=$contact_info{street}' -F 'city=$contact_info{city}' -F 'sub=$contact_info{state}' -F 'postal_code=$contact_info{zip}' -F 'country=$contact_info{country}' -F 'phone=$contact_info{phone}' -F 'fax=$contact_info{fax}' -F 'email=$paMeta->[15]' $sAuthorForm -F 'publication_status=unpublished' -F 'reference_title=$paMeta->[17]' -F 'same_reference_and_sequence_authors=yes' -F 'bioproject_id=$paMeta->[4]' -F 'biosample_id=$paMeta->[21]' -F 'submit=Create Template' $SUBMISSION_TEMPLATE > $sSbtFile";
	$sCmd =~ s/\r|\n//ig;
	printLogMsg($DEBUG, "Preparing to run - \n$sCmd\n");
	system($sCmd);

	if(-e $sSbtFile) {
		open($fhRead, "< $sSbtFile") or printLogMsg($ERROR, "Could not open file $sSbtFile for reading. Reason : $!");
		while(<$fhRead>) {
			chomp($_);
			if($_ =~ /^Submit-block/) {
				printLogMsg($DEBUG, "INFO : $sSbtFile file created successfully");
			} else {
				printLogMsg($ERROR, "ERROR : $sSbtFile file was created errorneously.");
			}
			last;
		}	
	} else {
		printLogMsg($ERROR, "ERROR : $sSbtFile file creation failed. Reason : $?");
	}
}

####################################################################################################################################################
# Description   : Used to create the file path on NCBI SRA FTP server using the SRX/SRR/SRS/SRP id  
# Parameters    : sRunId
#		  sRunId - SRA compatible experiment id or run id or sample id or study id for the read files to be downloaded
# Returns       : sFile
#		  sFile - NCBI SRA FTP file path for the specified id
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
# Description   : Used to create the file path on NCBI SRA FTP server using the SRX/SRR/SRS/SRP id  
# Parameters    : sRunId
#		  sRunId - SRA compatible experiment id or run id or sample id or study id for the read files to be downloaded
# Returns       : sFile
#		  sFile - NCBI SRA FTP file path for the specified id
# Modifications :

sub prepareCmt {
	my ($paMeta, $sOutDir) = @_;
	my $sSubName = (caller(0))[3];
	my $sCmtFile;
	my $fhWrite;

	$sCmtFile = $sOutDir."/".$paMeta->[0].".cmt";
	open($fhWrite, "> $sCmtFile") or printLogMsg($ERROR, "Could not open file $sCmtFile for writing. Reason : $!");
	
	print $fhWrite "StructuredCommentPrefix\t##Genome-Assembly-Data-START##\n";
	print $fhWrite "Assembly Method\t$paMeta->[11]\n";
	print $fhWrite "Genome Coverage\t$paMeta->[12]\n";
	print $fhWrite "Sequencing Technology\t$paMeta->[13]\n";
	print $fhWrite "StructuredCommentSuffix\t##Genome-Assembly-Data-END##";
	close($fhWrite);

	if((-e $sCmtFile) && (-s $sCmtFile)) {
		printLogMsg($DEBUG, "INFO : $sCmtFile file created successfully");
	} else {
		printLogMsg($ERROR, "ERROR : $sCmtFile file creation failed. Reason : $?");
	}
}


####################################################################################################################################################
# Description   : Used to download read files from NCBI SRA FTP server using wget
# Parameters    : phCmdLineArgs 
#		  phCmdLineArgs - reference to hash of command line arguments passed to the perl script
# Returns       : NA
# Modifications :

sub readInput {
	my ($sFile, $paMeta) = @_;
	my $nI;
	my $sLine;
	my $fhRead;
    	my $sSubName = (caller(0))[3];

	open($fhRead, "< $sFile") or printLogMsg($ERROR, "ERROR : Could not open $sFile file for reading. Reason : $!");
	while($sLine=<$fhRead>) {
		chomp($sLine);
		next if($sLine =~ /^#/);
		next if($sLine =~ /^\s+$/);
		@{$paMeta} = split(/\t/, $sLine);
		# @meta : This script needs columns 0, 5-20, and either one (or both) of 4 and 21
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
#		[18]= List of deprecated/bad EC-numbers
#		[19]= Isolation source of bacterial sample
#		[20]= Path to contact person address info
#		[21]= Biosample ID

		if(!length($paMeta->[0])) {
			printLogMsg($ERROR, "ERROR : Database name missing for $sLine.");
		}

		for($nI = 5; $nI <= 19; $nI++) {
			if(!length($paMeta->[$nI])) {
				printLogMsg($WARN, "WARNING : Missing meta data value for column ".($nI + 1)." field. Thus resulting files .sbt and .cmt would have missing information");
			}	
		}
		if (!length($paMeta->[21]) && !length($paMeta->[4])) {
			printLogMsg($WARN, "WARNING : Neither Biosample ID nor Bioproject ID detected.  Highly recommended to have at least one for the .cmt and .sbt files");
		}

		if (!length($paMeta->[20]) ) {
			printLogMsg($WARN, "WARNING : Path for author contact info not supplied.  Will use default parameters");
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

	if(exists($phCmdLineArgs->{'log'})) {
		open($logfh, "> $phCmdLineArgs->{'log'}") or die "Could not open $phCmdLineArgs->{'log'} file for writing.Reason : $!\n";
	}
	@aRequired = qw(input_file output_dir);
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

download_sra.pl - Script to download read files from NCBI SRA FTP site

=head1 SYNOPSIS

# USAGE : perl download_sra.pl -r <SRA id> -f <FTP server> -o <path to output dir> [ -u <FTP server username> -p <FTP server password> -n <number of retries> -l <path to log file> ]

	parameters in [] are optional

=head1 OPTIONS

	-r <run_id>	:	NCBI SRA compatible 9-character id. Could be Study id (SRPXXXXXX), Experiment id (SRXXXXXXX), Run id (SRRXXXXXX) or Sample id (SRSXXXXXX). X stands for digit. Mandatory
	
	-f <ftp>	:	Name of the NCBI FTP server. Currently it is ftp-trace.ncbi.nih.gov. Mandatory

	-o <output_dir>	:	Path to the output directory where the files will be downloaded. Mandatory

[	
	-u <username>	:	Username for the FTP server. Default: anonymous. Optional

	-p <password>	:	Password for the FTP server. Default: anonymous. Optional

	-n <num_retry>	:	Number of retries to download a file. wget default is 20. Optional

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
