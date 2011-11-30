#!/usr/local/bin/perl -w

#########################################################################################
#											#
# Name	      :	extract_snp_regions.pl							#
# Version     :	1.0									#
# Project     :	Prokaryotic SNP Verification Pipeline					#
# Description : Script to extract reference sequence surrounding predicted SNP position #
# Author      : Sonia Agrawal								#
# Date        :	November 4, 2011							#
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

#&parse_ref_genbank($cmdLineArgs{'ref_genbank'}, $cmdLineArgs{'output_dir'});

$snp_pos = &parse_snp_positions($cmdLineArgs{'snp_positions'});

&extract_seq($cmdLineArgs{'ref_genbank'}, $snp_pos, $flanking_bases, $cmdLineArgs{'output_dir'});

###############
# SUBROUTINES #
###############

# Description   : 
# Parameters    : Reference genome's genBank file path
#		: Output directory to store reference genome annotation file
# Returns       : NA
# Modifications :

sub parse_ref_genbank {
	my ($genbank, $result_dir) = @_;
	my ($ref_base,$ref_dir,$ref_ext) = fileparse($genbank,qr/\.[^.]*/);
	my $annot_file = "$result_dir/".$ref_base.".coords";
	open(FTAB, "> $annot_file") or die "Could not open $annot_file for writing.\nReason : $!\n";
	my $seqio_obj = Bio::SeqIO->new(-file => "$genbank", -format => "GenBank" );
	while (my $seqobj = $seqio_obj->next_seq()) {
		for my $feat_object ($seqobj->get_SeqFeatures) {
			if ($feat_object->primary_tag eq "CDS") {
				for my $locus ($feat_object->get_tag_values("locus_tag")) {
					$locus =~ s/\_//;
					print FTAB "$locus\t";
				}
				if ($feat_object->location->strand == -1) {
					print FTAB $feat_object->location->end."\t".$feat_object->location->start."\t";
				} else {
					print FTAB $feat_object->location->start."\t".$feat_object->location->end."\t";
				}
				for my $product ($feat_object->get_tag_values("product")) {
					print FTAB "$product\n";
				}
			}
		}
	}
	close(FTAB);
#	while(<FR>) {
#		my $line = $_;
#		chomp($line);
#		next if ($line =~ /^\s*$/);
#		next if ($line =~ /^#/);
#		my ($ref_id, $ref_path) = split(/\s+/,$line,2);
#		if (-e $ref_path) {
#			$ref_data->{$ref_id} = $ref_path;
#		} else {
#			die "File $ref_path does not exist. Check the reference file paths in $cmdLineArgs{'ref_info'}\n";
#		}
#	}
#	close(FR);
#	return $ref_data;
}

####################################################################################################################################################

# Description   :
# Parameters    : Path to the SNP positions file
# Returns       : Hash reference having snp position in reference genome as keys and rest of the 
#		  information in the snp_positions file as values	
# Modifications :

sub parse_snp_positions {
	my ($snp_file) = @_;
	my $snp_data;
	open(FR,"< $snp_file") or die "Could not open $snp_file file for reading.\nReason : $!\n";
	readline(FR);
	while(<FR>) {
		my $line = $_;
		chomp($line);
		next if ($line =~ /^\s*$/);
		next if ($line =~ /^#/);
		my ($refpos, $rest) = split(/\s+/,$line,2);
		$snp_data->{$refpos} = $rest;
	} 
	return $snp_data;
}

####################################################################################################################################################

# Description   :
# Parameters    : Reference genome's GenBank file path
#		: Hash of SNP positions
#		: Number of flanking bases required 
#		: Output directory to store file containing extracted SNP regions
# Returns       : NA	
# Modifications :

sub extract_seq {
	my ($genbank, $snp_data, $flank_num, $result_dir) = @_;
	my $seqio_obj = Bio::SeqIO->new(-file => "$genbank", -format => "GenBank" );
	my $seq_obj = $seqio_obj->next_seq;
	my $ref_seq = $seq_obj->seq();
	my $ref_length = length($ref_seq);
	if($ref_length <= 0) {
		die "GenBank file $genbank does not conatin the reference FASTA sequence\n";
	}
	my ($ref_base,$ref_dir,$ref_ext) = fileparse($genbank,qr/\.[^.]*/);
	my $seq_file = "$result_dir/extracted_snps_".$ref_base.".fna";
	open(FW, "> $seq_file") or die "Could not open $seq_file for writing.\nReason : $!\n";
	foreach my $pos (sort { $a <=> $b } keys %$snp_data) {
		print FW ">SNP_".$pos."\n";
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
			die "ERROR! : Required option $option not passed\n";
		}
	}
	if(exists($cmdLineArgs{'flanking_bases'})) {
		$flanking_bases = $cmdLineArgs{'flanking_bases'};	
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
