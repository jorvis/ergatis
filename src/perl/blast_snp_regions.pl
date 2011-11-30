#!/usr/local/bin/perl -w

#########################################################################################
#											#
# Name	      :	blast_snp_regions.pl							#
# Version     :	1.0									#
# Project     :	Prokaryotic SNP Verification Pipeline					#
# Description : Script to BLAST reference sequence surrounding predicted SNP position	#
#		against query genome							#
# Author      : Sonia Agrawal								#
# Date        :	November 9, 2011							#
#											#
#########################################################################################

use strict;
# Extended processing of command line options
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
# Prints usage message from embedded pod documentation
use Pod::Usage;
# Include packages here 
use Bio::SeqIO;
use File::Basename;

#############
# CONSTANTS #
#############


###########
# GLOBALS #
###########
my %cmdLineArgs = ();
my $evalue = 0.0001;
my $filter = "F";
my $blast_prog = "blastn";
# Log file handle;
my $logfh;
my ($ERROR, $WARN, $DEBUG) = (1, 2, 3);

################
# MAIN PROGRAM #
################
GetOptions(\%cmdLineArgs,
	   'blast_exec|b=s',
	   'database|d=s',
	   'query_list|q=s',
	   'evalue|e=s',
	   'filter|f=s',
	   'other_args|a=s',
	   'output_dir|o=s',
	   'log|l=s', 
	   'help|h'
	  ) or pod2usage(2);

pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDERR} ) if ($cmdLineArgs{'help'});

&checkCmdLineArgs();

open(FR,"< $cmdLineArgs{'query_list'}") or die "Could not open file $cmdLineArgs{'query_list'} for reading.\nReason : $!\n";
while(<FR>) {
	my $file = $_;
	chomp($file);
	next if ($file =~ /^\s*$/);
	if(-e $file) {
		my ($file_base,$file_dir,$file_ext) = fileparse($file,qr/\.[^.]*/);
		$file_base =~ /extracted_snps_(.+)/;
		my $ref_name = $1;
		&blast_snps($file,$ref_name);
	}
}
close(FR);

###############
# SUBROUTINES #
###############

# Description   : 
# Parameters    : Reference genome's genBank file path
#		: Output directory to store reference genome annotation file
# Returns       : NA
# Modifications :

sub blast_snps {
	my ($blast_query,$chr) = @_;
	my ($db_base,$db_dir,$db_ext) = fileparse($cmdLineArgs{'database'},qr/\.[^.]*/);
	my $blast_db = $db_dir.$db_base;
	my $blast_output = $cmdLineArgs{'output_dir'}."/".$db_base."_refmol_".$chr.".raw";
	my $blast_command = "$cmdLineArgs{'blast_exec'} -p $blast_prog -F $filter -e $evalue -d $blast_db -i $blast_query -o $blast_output $cmdLineArgs{'other_args'}";
	print STDOUT "BLAST COMMAND : $blast_command\n";
	if(system($blast_command)) {
		 die "Could not execute the BLAST command $blast_command\n";
	}
}

####################################################################################################################################################

# Description   : Used to check the correctness of the command line arguments passed to the script. The script exits if required arguments are missing. 
# Parameters    : NA
# Returns       : NA
# Modifications :

sub checkCmdLineArgs {
	my @required = qw(blast_exec database query_list output_dir);
	foreach my $option(@required) {
		if(!defined($cmdLineArgs{$option})) {
			die "ERROR! : Required option $option not passed\n";
		}
	}
	if(exists($cmdLineArgs{'evalue'})) {
		$evalue = $cmdLineArgs{'evalue'};	
	}
	if(exists($cmdLineArgs{'filter'})) {
		$filter = $cmdLineArgs{'filter'};
	}
	if(exists($cmdLineArgs{'log'})) {
		open($logfh, "> $cmdLineArgs{'log'}") or die "Could not open $cmdLineArgs{'log'} file for writing.Reason : $!\n"
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
		die "$msg\n" if($level == $ERROR);
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
