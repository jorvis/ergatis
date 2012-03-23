#!/usr/local/bin/perl -w

#################################################################################################
#												#
# Name	      	: extract_snp_regions.pl							#
# Version     	: 2.0										#
# Project     	: Prokaryotic SNP Verification Pipeline						#
# Description 	: Script to extract reference sequence surrounding predicted SNP position 	#
# Modifications	: 2.0 - Accomodate for change in SNP panel format				#
# Author      	: Sonia Agrawal									#
# Date        	: November 4, 2011								#
#		  January 11, 2012		  						#
#												#
#################################################################################################

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
my $flanking_bases = 20;
my $snp_pos;
# Log file handle;
my $logfh;
my ($ERROR, $WARN, $DEBUG) = (1, 2, 3);

################
# MAIN PROGRAM #
################
GetOptions(\%cmdLineArgs,
	   'ref_genbank|r=s',
	   'snp_positions|s=s',
	   'flanking_bases|f=i',
	   'output_dir|o=s',
	   'log|l=s', 
	   'help|h'
	  ) or pod2usage(2);

pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDERR} ) if ($cmdLineArgs{'help'});

&checkCmdLineArgs();

$snp_pos = &parse_snp_positions($cmdLineArgs{'snp_positions'});

&extract_seq($cmdLineArgs{'ref_genbank'}, $snp_pos, $flanking_bases, $cmdLineArgs{'output_dir'});

###############
# SUBROUTINES #
###############

# Description   : Parse SNP panel files to create a unique SNP panel
# Parameters    : snp_file = Path to the SNP positions file
# Returns       : snp_data = Hash reference having snp position in reference genome as keys and rest of the 
#		  information in the snp_positions file as values	
# Modifications : The return hash reference has molecule (LOCUS id) as first key and then snp position as second key
#		  Accepts path to a list of SNP panels and creates a merged hash of all the SNPs in all the panels in the list

sub parse_snp_positions {
	my ($snp_file) = @_;
	my $snp_data;
	open(SL,"< $snp_file") or printLogMsg($ERROR, "ERROR! : Could not open $snp_file file for reading.Reason : $!");
	# Loop through SNP panel files list
	while(my $file = <SL>) {
		chomp($file);
		next if($file =~ /^\s*$/);
		next if($file =~ /^#/);
		# Foreach SNP panel file extract the SNP position and molecule name
		if(-e $file) {
			open(FR,"< $file") or printLogMsg($ERROR, "ERROR! : Could not open $file file for reading.Reason : $!");
			readline(FR);
			while(<FR>) {
				my $line = $_;
				chomp($line);
				next if ($line =~ /^\s*$/);
				next if ($line =~ /^#/);
				my ($mol, $refpos, $rest) = split(/\s+/, $line, 3);
				$snp_data->{$mol}{$refpos} = $rest;
			}	
			close(FR); 
		}
	}
	close(SL);
	return $snp_data;
}

####################################################################################################################################################

# Description   : Extracts sub-sequence from the reference genome FASTA sequence SNP position along with flanking bases upstream and downstream of it
# Parameters    : genbank = Reference genome's GenBank file path
#		: snp_data = Hash of SNP positions
#		: flank_num = Number of flanking bases required 
#		: result_dir = Output directory to store file containing extracted SNP regions
# Returns       : NA	
# Modifications : Added code to accomodate for multiple reference GenBank files and extracting sequences from them based on information from 
#		  SNP panel file ie. LOCUS id of the SNP position

sub extract_seq {
	my ($genbank, $snp_data, $flank_num, $result_dir) = @_;
#	open(RL, "< $genbank") or printLogMsg(1,"ERROR! : Could not open $genbank file for reading.Reason : $!");
	my $seq_file = "$result_dir/extracted_snps.fna";
	open(FW, "> $seq_file") or printLogMsg(1,"ERROR! : Could not open $seq_file for writing.\nReason : $!");
#	while(my $file = <RL>) {
#		chomp($file);
#		next if($file =~ /^\s*$/);
#		next if($file =~ /^#/);
#		if(-e $file) {
			my $seqio_obj = Bio::SeqIO->new(-file => "$genbank", -format => "GenBank" );
			my $seq_obj = $seqio_obj->next_seq;
			my $locus = $seq_obj->display_id();
			if(exists($snp_data->{$locus})) {
				my $ref_seq = $seq_obj->seq();
				my $ref_length = length($ref_seq);
				if($ref_length <= 0) {
					printLogMsg($ERROR,"ERROR! : GenBank file $genbank does not conatin the reference FASTA sequence. Unable to extract sequences from $genbank");
				} else {
					foreach my $pos (sort { $a <=> $b } keys %{$snp_data->{$locus}}) {
						print FW ">".$locus."_SNP_".$pos."\n";
						my $left_limit = $pos - $flank_num;
						if ($left_limit < 1) {
							$left_limit = 1;
						}
						my $right_limit = $pos + $flank_num;
						if ($right_limit > $ref_length) {
							$right_limit = $ref_length;
						}
						my $snp_region = substr($ref_seq, ($left_limit - 1), ($right_limit - $left_limit + 1));
						print FW "$snp_region\n";
					}
				}
			} else {
				printLogMsg($ERROR, "ERROR! : SNPs for $locus is not found in the SNP panel. Unable to extract sequences for $locus from GenBank file $genbank");
			}
#		} else {
#			printLogMsg(2, "WARNING! : Reference GenBank file $file does not exist. Unable to extract sequences from $file");
#		}
#	}
	close(FW);
}

####################################################################################################################################################

# Description   : Used to check the correctness of the command line arguments passed to the script. The script exits if required arguments are missing. 
# Parameters    : NA
# Returns       : NA
# Modifications :

sub checkCmdLineArgs {
	my @required = qw(ref_genbank snp_positions output_dir);
	foreach my $option(@required) {
		if(!defined($cmdLineArgs{$option})) {
			printLogMsg(1,"ERROR! : Required option $option not passed");
		}
	}
	if(exists($cmdLineArgs{'flanking_bases'}) && defined($cmdLineArgs{'flanking_bases'})) {
		$flanking_bases = $cmdLineArgs{'flanking_bases'};	
	}
	if(exists($cmdLineArgs{'log'})) {
		open($logfh, "> $cmdLineArgs{'log'}") or printLogMsg(1,"Could not open $cmdLineArgs{'log'} file for writing.Reason : $!");
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
		die "\n" if($level == $ERROR);
	}	
}

__END__

#####################
# POD DOCUMENTATION #
#####################

=head1 NAME

extract_snp_regions.pl - Script to extract reference sequence surrounding predicted SNP positions

=head1 SYNOPSIS

USAGE : perl extract_snp_regions.pl --ref_genbank <reference genbank file> --snp_positions <snp_panel_file> --output_dir <output_dir> [ --flanking_bases <number of flanking bases> --log <log_file> ]

	parameters in [] are optional

=head1 OPTIONS

	-ref_genbank	: Path to reference GenBank file, for a chromosome. Mandatory   

	-snp_positions	: path to SNP panel file having reference positions as first column and LOCUS id from the reference GenBank file as the second column, may have additional columns. Mandatory

	-output_dir	: /path/to/output_dir. Mandatory

	-flanking_bases	: number of bases to be used as flanking bases on either side of the SNP base. Default is 20. Optional

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
