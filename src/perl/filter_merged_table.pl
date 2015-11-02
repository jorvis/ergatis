#!/usr/local/bin/perl -w

#################################################################################
#										#
# Name	      : filter_merged_table.pl						#
# Version     : 1.0								#
# Project     : SNP Verification pipeline					#
# Description : Script to remove "No Hits" and greater than 1 base results from #
#		merged table							#
# Author      : Sonia Agrawal							#
# Date        : April 4, 2012							#
#										#
#################################################################################

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
my %cmdLineArgs = ();
my %spg_df = ();
my %spg_ff = ();
my %head = ();
my %seq = ();
my @header = ();
my @rows = ();
my @discarded = ();
my @filtered = ();
my @merged_data = ();
my @genes = ();
my $query_list = [];
# Log file handle;
my $logfh;
my ($line, $resnum, $i, $no_hit_flag, $gene);
my ($filtered_file, $discarded_file, $fasta_file, $col_genename, $col_spg);
my $linenum = 1;
my ($ERROR, $WARN, $DEBUG) = (1, 2, 3);

################
# MAIN PROGRAM #
################
GetOptions(\%cmdLineArgs,
	   'merged_table|m=s',
	   'query_list|q=s',
	   'output_dir|o=s',
	   'log|l=s', 
	   'help|h'
	  ) or pod2usage();

pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDERR} ) if ($cmdLineArgs{'help'});

&checkCmdLineArgs();

$query_list = mapQueryList($cmdLineArgs{'query_list'});
$col_genename = 4 + @{$query_list};
$col_spg = $col_genename + 5;

# Reading input merged table
open(MT, "< $cmdLineArgs{'merged_table'}") or printLogMsg($ERROR, "ERROR! : Could not open $cmdLineArgs{'merged_table'} for reading");
@merged_data = <MT>;
close(MT);

# Output files
$filtered_file = $cmdLineArgs{'output_dir'}."/filtered_merged_table.txt";
$discarded_file = $cmdLineArgs{'output_dir'}."/discarded_merged_table.txt";
$fasta_file = $cmdLineArgs{'output_dir'}."/filtered_fasta_merged_table.fna";

foreach $line (@merged_data) {
	chomp $line;
	next if($line =~ /^\s*$/);
	@rows = split(/\t/,$line);
	if($line =~ /^molecule/) {
		@header = @rows;
		next;
	}
	$no_hit_flag = 0;
	for ($i=4; $i < $col_genename; $i++) {
		# Applying filter to remove rows form merged table where queries have "No Hit" or "indel" or more than one nucleotide at a SNP position
		if(($rows[$i] eq "No Hit") || (length($rows[$i]) > 1) || ($rows[$i] eq "indel")) {
			$no_hit_flag = 1;
			last;
		}
	}
	@genes = split(/\//,$rows[$col_genename]);
	if($no_hit_flag) {
		foreach $gene (@genes) {
			$spg_df{$rows[0]}{$gene} = 0 if(!exists($spg_df{$rows[0]}{$gene}));
			$spg_df{$rows[0]}{$gene}++;
		}
		# Pushing line numbers to be discarded from the merged table based on the above filter
		push(@discarded, $linenum);
	} else {
		foreach $gene (@genes) {
			$spg_ff{$rows[0]}{$gene} = 0 if(!exists($spg_ff{$rows[0]}{$gene}));
			$spg_ff{$rows[0]}{$gene}++;
		}
		push(@filtered, $linenum);
		# Storing sequence for each query genome in the filtered results to create FASTA file		
		for ($i=3; $i < $col_genename; $i++) {
			$seq{$header[$i]} .= $rows[$i];
		}
	}
	$linenum++;
}

# Printing filtered, discarded and FASTA files 
print_results($filtered_file, \@filtered, \%spg_ff);
print_results($discarded_file, \@discarded, \%spg_df);
print_fasta($fasta_file, \%seq);

###############
# SUBROUTINES #
###############

# Description   : Create a map file to map name of query genomes to the names of the contigs belonging to the query genome  
# Parameters    : query_list = Path to the list 
# Returns       : queries = Reference to an array of query genome names 
# Modifications :

sub mapQueryList {
	my ($query_file) = @_;
	my @queries = ();
	my ($file, $file_base, $file_dir, $file_ext);
	open(QL, "< $query_file") or printLogMsg($ERROR, "Could not open $query_file file for reading. Reason : $!");
	# Loop through query genome list file
	while($file = <QL>) {
		chomp($file);
		next if ($file =~ /^\s*$/);
		next if ($file =~ /^#/);
		if ((-e $file) && (-s $file > 0)) {
			($file_base,$file_dir,$file_ext) = fileparse($file,qr/\.[^.]*/);
			push @queries, $file_base;
		} else {
			printLogMsg($ERROR, "Specified $file file does not exist or is empty.");
		}
	}
	close(QL);
	return(\@queries);
}

####################################################################################################################################################

# Description   : Used to print FASTA sequence from the filtered results 
# Parameters    : filename - Path to the output file
#		  sequences - Reference to a hash containing query name and sequence of SNP bases to be printed
# Returns       : NA
# Modifications :

sub print_fasta {
	my ($filename, $sequences) = @_;
	my $s;
	open(OUT, "> $filename") or printLogMsg($ERROR, "ERROR! : Could not open $filename for writing");
	for($s=3; $s < $col_genename; $s++) {
		print OUT ">$header[$s] length:".length($sequences->{$header[$s]})."\n";
		print OUT "$sequences->{$header[$s]}\n";
	} 
	close(OUT);
}

####################################################################################################################################################

# Description   : Used to print results and adjust the snps_per_gene based on the lines_to_print array 
# Parameters    : filename - Path to the output file
#		  lines_to_print - Reference to an array containing the line numbers to be print from the merged table
#		  snps_count - Reference to a hash containing the gene name and the number of SNPs in that gene 
# Returns       : NA
# Modifications :

sub print_results {
	my ($filename, $lines_to_print, $snps_count) = @_;
	my @row = ();
	my @genes = ();
	my @snps_per_gene = ();
	my ($res, $gene);
	open(OUT, "> $filename") or printLogMsg($ERROR, "ERROR! : Could not open $filename for writing");
	chomp($merged_data[0]);
	print OUT "$merged_data[0]\n";
	foreach $res (@{$lines_to_print}) {
		chomp($merged_data[$res]);
		@row = split(/\t/,$merged_data[$res]);
		if($row[$col_genename] eq "intergenic") {
			$row[$col_spg] = "NA";
			print OUT join("\t", @row);
			print OUT "\n";
			next;
		}
		@snps_per_gene = ();
		@genes = split(/\//, $row[$col_genename]);
		# Updating snps_per_gene column based on the results to be printed
		foreach $gene (@genes) {
			if(exists($snps_count->{$row[0]}{$gene})) {
				push(@snps_per_gene,$snps_count->{$row[0]}{$gene});
			} else {
				printLogMsg($ERROR, "ERROR! : Could not find SNPs per gene for ".$row[0]." and $gene at position $row[1]");
			}
		}
		$row[$col_spg] = join("/", @snps_per_gene);
		print OUT join("\t", @row);
		print OUT "\n";		
	}
	close(OUT);
}

####################################################################################################################################################

# Description   : Used to check the correctness of the command line arguments passed to the script. The script exits if required arguments are missing. 
# Parameters    : NA
# Returns       : NA
# Modifications :

sub checkCmdLineArgs {
	my @required = qw(merged_table query_list output_dir);
	foreach my $option (@required) {
		if(!defined($cmdLineArgs{$option})) {
			printLogMsg($ERROR,"ERROR! : Required option $option not passed");
		}
	}
	if(exists($cmdLineArgs{'log'})) {
		open($logfh, "> $cmdLineArgs{'log'}") or printLogMsg($ERROR, "ERROR! : Could not open $cmdLineArgs{'log'} file for writing.Reason : $!");
	}
}

####################################################################################################################################################

# Description   : Used to handle logging of messages(errors and warnings) during the execution of the script
# Parameters    : level = can be ERROR, WARNING or INFO
#		  msg   = msg to be printed in the log file or to STDERR
# Returns       : NA
# Modifications : 22

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

filter_merged_table.pl - Script to filter the merged table obtained after SNP verification to remove rows with No Hits or indels or more than one nucleotide at a SNP position in a query genome.
			 It also prints the FASTA file of all the SNPs from the filtered results

=head1 SYNOPSIS

USAGE : filter_merged_table.pl -m <merged_table> -q <num_of_query> -o <output_directory> [ -l <log_file> ]

	parameters in [] are optional

=head1 OPTIONS

	-m <merged_table>	: Path to merged table obtained after SNP verification. Mandatory

	-q <query_num>		: Number of query genomes in the merged table. Mandatory

	-o <output_dir>		: Path to output directory where all the result files will be printed. Mandatory

	-l <log_file>		: Path to log file. Optional	

=head1 DESCRIPTION



=head1 INPUT

	Refer to snp_verify to see the required format of merged table. If the exact format is not present, the script may print incorrect results 

=head1 OUTPUT

	filtered_merged_table.txt - Filtered merged table without rows containing "No Hits" or "indel" or more than one nucleotide at a SNP position in a query genome

	discarded_merged_table.txt - All the rows removed from merged table to obtain the above filtered_merged_table.txt

	filtered_fasta_merged_table.fna - Multi-FASTA file containing the sequence of SNP nucleotides for each query genome in the merged table.	

=head1 AUTHOR

	Sonia Agrawal
	Bioinformatics Software Engineer II
	Institute for Genome Sciences
	School of Medicine
	University of Maryland
	Baltimore, MD - 21201
	sagrawal@som.umaryland.edu

==cut
