#!/usr/local/bin/perl -w

#########################################################################
#									#
# Name	      :	configure_manatee.pl					#
# Version     :	1.0							#
# Project     :	CloVR-microbe pipeline					#
# Description : Script to create MySQL user and database to store	# 
#		prok annotation during the pipeline run			#
# Author      : Sonia Agrawal						#
# Date        :	April 15, 2015						#
#									#
#########################################################################

use strict;
# Extended processing of command line options
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
# Prints usage message from embedded pod documentation
use Pod::Usage;
# Include packages here 
use DBI;

###########
# GLOBALS #
###########
my %cmdLineArgs;
my $mysql_username = 'root';
#my $mysql_password = 'Cl0vR_mY$?L';
my $mysql_password = 'clovr_mysql';
my $host = 'localhost';
# MySQL connection handle
my $dbh;
# Log file handle;
my $logfh;
my ($ERROR, $WARN, $DEBUG) = (1, 2, 3);

################
# MAIN PROGRAM #
################
GetOptions(\%cmdLineArgs,
	   'db_name|d=s',
	   'db_username|u=s',
	   'db_password|p=s',
	   'db_server|s=s',
	   'action|a=s',
	   'log|l=s', 
	   'debug|n=i',
	   'help|h'
	  ) or pod2usage(2);

pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDERR} ) if ($cmdLineArgs{'help'});

&checkCmdLineArgs();

# Connect to MySQL server as root
$dbh = DBI->connect("DBI:mysql:host=$host",$mysql_username,$mysql_password) or printLogMsg($ERROR,"Could not connect to $host as $mysql_username");
printLogMsg($DEBUG,"Connection to MySQL successful");

if($cmdLineArgs{'action'} eq "create") {
# Create database
	$dbh->do("CREATE DATABASE IF NOT EXISTS $cmdLineArgs{'db_name'}") or printLogMsg($ERROR,"Could not create $cmdLineArgs{'db_name'} on $host");
	printLogMsg($DEBUG,"Creation of $cmdLineArgs{'db_name'} successful");

# Create DB username 
#	$dbh->do("CREATE USER '$cmdLineArgs{db_username}'\@'localhost' IDENTIFIED BY '$cmdLineArgs{db_password}'") or printLogMsg($ERROR,"Could not create MySQL user $cmdLineArgs{'db_username'}");
#	printLogMsg($DEBUG,"Creation of $cmdLineArgs{'db_username'} successful");


# Grant permissions on the database to the user
# User will be automatically created if they do not currently exist
	$dbh->do("GRANT ALL PRIVILEGES ON ".$cmdLineArgs{'db_name'}.".* TO '$cmdLineArgs{db_username}'\@'localhost' IDENTIFIED BY '$cmdLineArgs{db_password}'") or printLogMsg($ERROR,"Could not grant permissions to $cmdLineArgs{'db_username'} on $cmdLineArgs{'db_name'}");
	printLogMsg($DEBUG,"Grant permissions to $cmdLineArgs{'db_username'} on $cmdLineArgs{'db_name'} successful");

} elsif($cmdLineArgs{'action'} eq "drop") {
# Drop database	
	$dbh->do("DROP DATABASE $cmdLineArgs{'db_name'}") or printLogMsg($ERROR,"Could not drop $cmdLineArgs{'db_name'} on $host");
	printLogMsg($DEBUG,"Deleting of $cmdLineArgs{'db_name'} successful");
# Drop user
	$dbh->do("DROP USER '$cmdLineArgs{db_username}'\@'localhost'") or printLogMsg($ERROR,"Could not drop $cmdLineArgs{'db_username'} on $host");
	printLogMsg($DEBUG,"Deleting of $cmdLineArgs{'db_username'} successful");
}

###############
# SUBROUTINES #
###############

# Description   : Used to check the correctness of the command line arguments passed to the script. The script exits if required arguments are missing. 
# Parameters    : NA
# Returns       : NA
# Modifications :

sub checkCmdLineArgs {
	if(exists($cmdLineArgs{'log'})) {
		open($logfh, "> $cmdLineArgs{'log'}") or printLogMsg($ERROR,"Could not open $cmdLineArgs{'log'} file for writing.Reason : $!");
	}
#	my @required = qw(db_name db_username db_password action);
	my @required = qw(db_name db_username db_password action);
	for my $option(@required) {
		if(!defined($cmdLineArgs{$option})) {
			printLogMsg($ERROR,"ERROR! : Required option $option not passed");
		}
	}
	if(defined($cmdLineArgs{'db_server'})) {
		$host = $cmdLineArgs{'db_server'};
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

mysql_credential.pl - Script used to create MySQL user and database to store prokaryotic annotation durind the pipeline run

=head1 SYNOPSIS

USAGE :	mysql_credential.pl --c <flag> --u <DB username> --p <DB password> --d <DB name> [--l <log file>]

	parameters in [] are optional

=head1 OPTIONS

	--c <flag>		= 0 if DB credentials do not exist and need to create username and password
			  	  1 if DB credentials exist 
	
	--u <DB username>	= Existing or desired DB username

	--p <DB passoword>	= Existing or desired DB password

	--d <DB name>		= Name of the new DB to store the genome annotation

	--l <log file>		= /path/to/log_file. Optional

=head1 DESCRIPTION

	The program creates a new database to store genome annotation produced during the prokaryotic pipeline execution.
	It also creates the database user and password if --c parameter is set to 0. It grants appropriate permissions on 
	the new DB to the DB user. 

=head1 AUTHOR

	Sonia Agrawal
	Bioinformatics Software Engineer II
	Institute for Genome Sciences
	School of Medicine
	University of Maryland
	Baltimore, MD - 21201
	sagrawal@som.umaryland.edu

==cut
