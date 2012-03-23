#!/usr/local/bin/perl -w

#########################################################################################
#											#
# Name	      :	snp_verify.pl								#
# Version     : 2.0									#
# Project     : SNP Verification							#
# Description : Script to parse BLAST results and annotate verified SNPs.		#
# Author      : Sonia Agrawal								#
# Date        : February 6, 2012							#
#											#
#########################################################################################

use strict;
# Extended processing of command line options
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
# Prints usage message from embedded pod documentation
use Pod::Usage;
# Include packages here 
# Creates a tree for searching SNP positions on genes
use IntervalTree;
# Parses GenBank file 
use Bio::SeqIO;
# Parses BLAST results
use Bio::SearchIO;
use File::Basename;

###########
# GLOBALS #
###########
my %cmdLineArgs = ();
my %snps_per_gene = ();
my $snp_panel = {};
my $query_map = {};
my $results = {};
my $coords = {};
my $gbk_seq = {};
my $query_list = [];
my ($head_labels, $molecule, $refpos, $snp_tag, $refbase, $ol, $gene_id, $molname, $gene_name, $pos_in_gene, $geneseq, $pos_in_codon, $refcodon, $refaa, $qaa, $qcodon, $syn, $q, $num_hits, $snp_string, $num);
my (@gene_annot, @ols);
my ($output_file, $temp_file, $table);
my $flag = 0;
# Object for building tree of gene start and stop
my $iTree = new IntervalTree;
# Log file handle;
my $logfh;
my ($ERROR, $WARN, $DEBUG) = (1, 2, 3);

################
# MAIN PROGRAM #
################
GetOptions(\%cmdLineArgs,
	   'raw_blast_list|b=s',
	   'ref_genbank|r=s',
	   'snp_positions|s=s',
	   'query_list|q=s',
	   'flanking_bases|f=i',
	   'output_dir|o=s',
	   'log|l=s', 
	   'help|h'
	  ) or pod2usage();

pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDERR} ) if ($cmdLineArgs{'help'});

checkCmdLineArgs();

printLogMsg($DEBUG, "INFO :: Parsing SNP panels list file to find a unique list of SNP positions");
$snp_panel = parseSnpPanels($cmdLineArgs{'snp_positions'});

printLogMsg($DEBUG, "INFO :: Parsing queries list file to create the mapping file of query contigs to queryname");
($query_map, $query_list) = mapQueryList($cmdLineArgs{'query_list'});

printLogMsg($DEBUG, "INFO :: Parsing raw BLAST results list file");
$results = parseRawBlast($cmdLineArgs{'raw_blast_list'}, $cmdLineArgs{'flanking_bases'}, $query_map);

printLogMsg($DEBUG, "INFO :: Parsing reference genBank list file");
($coords, $gbk_seq) = parseGenbankFiles($cmdLineArgs{'ref_genbank'});

printLogMsg($DEBUG, "INFO :: Building tree to search for location of SNPs within genes");
$iTree->buildTree();

$output_file = $cmdLineArgs{'output_dir'}."/merged_table.txt";
# Temporary intermediate output file created to store SNP and its annotation 
$temp_file = $cmdLineArgs{'output_dir'}."/merged_table.temp";
open(MT, "> $temp_file") or printLogMsg($ERROR, "ERROR! : Could not open $temp_file file for writing.Reason : $!");

# Creating output file header
$head_labels = makeHeader($query_list);
print MT $head_labels."\n";

# Looping through the SNP panel hash with chromosome as first key and reference SNP position as second key
foreach $molecule (keys %$snp_panel) {
	foreach $refpos (sort {$a <=> $b} keys %{$snp_panel->{$molecule}}) {
		# Initialize gene annotation array for each SNP position. 
		#(syn/nsyn, gene_name,product,gene_start,gene_stop,gene_length,snps_per_gene,pos_in_gene,ref_codon,ref_aa,query_codon,query_aa)  
		@gene_annot = ('','','','','','','','','','','','');
		$flag = 0;
		$snp_tag = $molecule."_SNP_".$refpos;
		print MT "$molecule\t$refpos\t";	
		# Fetching reference base from the reference GenBank file for the reference SNP position
		$refbase = getSeqFromParent($gbk_seq, $molecule, $refpos, $refpos, 1);
		if(length($refbase) == 0) {
			printLogMsg($WARN, "WARNING! : No base found in $molecule at $refpos, skipping $refpos in the output file");
			next;
		}
		# Searching the tree to find if the reference SNP position is within gene(s)
		@ols = $iTree->searchInterval( ($refpos - 1), $refpos );
		# If reference SNP position is with gene(s)
		if(@ols > 0) {
			
			# If overlapping genes then more than one gene name would be returned from the tree search, so looping through the search results
			# Results for more than one gene per SNP will be printed as "/" separated 
			foreach $ol (@ols) {
				$gene_id = $ol->[2];
				($molname,$gene_name) = split(/\|/,$gene_id,2);
				# Counter hash for calculating SNPs per gene
				$snps_per_gene{$molecule}{$gene_name} = 0 if(!exists($snps_per_gene{$molecule}{$gene_name}));
				$snps_per_gene{$molecule}{$gene_name}++;
				# Next if gene returned belongs to another chromosome than the one being searched.
				next if ($molname ne $molecule);
				$flag = 1;
				# Gene name - locus tag from reference GenBank file of the genome 
				$gene_annot[1] .= ($gene_annot[1] eq "")? $gene_name : "/".$gene_name;
				# Function of the gene - product tag from reference GenBank file of the genome
				$gene_annot[2] .= ($gene_annot[2] eq "")? $coords->{$molname}{$gene_name}{'product'} : "/".$coords->{$molname}{$gene_name}{'product'};
				# Start position of the gene within the reference genome
				$gene_annot[3] .= ($gene_annot[3] eq "")? $coords->{$molname}{$gene_name}{'start'} : "/".$coords->{$molname}{$gene_name}{'start'};
				# End position of the gene within the reference genome
				$gene_annot[4] .= ($gene_annot[4] eq "")? $coords->{$molname}{$gene_name}{'end'} : "/".$coords->{$molname}{$gene_name}{'end'};
				# Length of the gene
				$gene_annot[5] .= ($gene_annot[5] eq "")? $coords->{$molname}{$gene_name}{'length'} : "/".$coords->{$molname}{$gene_name}{'length'};
				# Calculate relative position of the SNP in the gene based on strand information of the reference gene
				if($coords->{$molname}{$gene_name}{'strand'} == 1) {
					$pos_in_gene = $refpos - $coords->{$molname}{$gene_name}{'start'} + 1;
					# Sequence of the gene 5' - 3'
					$geneseq = getSeqFromParent($gbk_seq, $molecule, $coords->{$molname}{$gene_name}{'start'}, $coords->{$molname}{$gene_name}{'end'}, $coords->{$molname}{$gene_name}{'strand'});
					$gene_annot[7] .= ($gene_annot[7] eq "")? $pos_in_gene : "/".$pos_in_gene ;
				} elsif($coords->{$molname}{$gene_name}{'strand'} == -1) {
					$pos_in_gene = $coords->{$molname}{$gene_name}{'start'} - $refpos + 1;
					# Sequence of the gene 3' - 5'
					$geneseq = getSeqFromParent($gbk_seq, $molecule, $coords->{$molname}{$gene_name}{'end'}, $coords->{$molname}{$gene_name}{'start'}, $coords->{$molname}{$gene_name}{'strand'});
					$gene_annot[7] .= ($gene_annot[7] eq "")? $pos_in_gene : "/".$pos_in_gene ;
				} else {
					printLogMsg($ERROR, "ERROR ! : Strand information is incorrect for $molname-$gene_name");
				}
				# Relative position of the reference SNP in the codon. It can be 1,2 or 0
				$pos_in_codon = $pos_in_gene % 3;
				# 0 means 3rd base in the codon
				$pos_in_codon = 3 if( $pos_in_codon == 0);
				# Calculate reference codon using relative position in gene and in codon
				$refcodon = substr($geneseq, ($pos_in_gene - $pos_in_codon), 3);
				$gene_annot[8] .= ($gene_annot[8] eq "")? $refcodon : "/".$refcodon;
				# Calculate reference amino acid
				$refaa = translateCodon($refcodon);
				$gene_annot[9] .= ($gene_annot[9] eq "")? $refaa : "/".$refaa;
				# Find query codon by replacing query base in the reference codon. Compute query amino acid. 
				# Check if this position has synonymous or non-synonymous SNP by comparing reference amino acid and query amino acid
				($qcodon, $qaa, $syn) = getSnpType($pos_in_codon, $refcodon,$refbase,$coords->{$molname}{$gene_name}{'strand'}, $snp_tag, $results);
				$gene_annot[0] .= ($gene_annot[0] eq "")? $syn : "/".$syn;
				$gene_annot[10] .= ($gene_annot[10] eq "")? $qcodon : "/".$qcodon;
				$gene_annot[11] .= ($gene_annot[11] eq "")? $qaa : "/".$qaa;
			}
		} 
		if (!$flag) {
			@gene_annot = ('NA','intergenic','intergenic','NA','NA','NA','NA','NA','NA','NA','NA','NA');
		}
		# Print type of SNP (syn or non-syn) and reference base
		print MT "$gene_annot[0]\t";
		print MT "$refbase\t";	
		# Calculate number of hits in each query and print query base
		$num_hits = "";
		# If no hit found for that SNP pposition in any of the queries
		if(!exists($results->{$snp_tag})) {
			$num_hits = "";
			map {
				print MT "No Hit\t";
				$num_hits .= "0\t";
			} @{$query_list};
			printLogMsg($DEBUG, "INFO: No query blast hits at SNP location $snp_tag");	
		} else {
			foreach $q ( @{$query_list} ) {
				$snp_string = "No Hit";
				if(exists($results->{$snp_tag}{$q})) {
					$num = 0;
					$snp_string = join("", keys %{$results->{$snp_tag}{$q}});
					$num += values %{$results->{$snp_tag}{$q}};
					$num_hits .= "$num\t";
				} else {
					$num_hits .= "0\t";
				}
				print MT "$snp_string\t";
			}
		}
		# Print all the annotation for the gene and number of hits.
		print MT "$gene_annot[1]\t$gene_annot[2]\t$gene_annot[3]\t$gene_annot[4]\t$gene_annot[5]\t\t$gene_annot[7]\t$gene_annot[8]\t$gene_annot[9]\t$gene_annot[10]\t$gene_annot[11]\t$num_hits"."$snp_panel->{$molecule}{$refpos}\n";
	}
}
close(MT);

# Add snps_per_gene column to the final output table
# First printed everything else in a temp file and then added this column
if(-e $temp_file) {
	my $col_spg = 4 + @{$query_list} + 5;
	print STDERR "colspg = $col_spg\n";
	printLogMsg($DEBUG, "INFO :: Adding SNPs per gene information to the final output file");
	printLogMsg($DEBUG, "INFO :: Printing output file $output_file");
	addSnpsPerGene($col_spg, $temp_file, $output_file, \%snps_per_gene);
} else {
	printLogMsg($ERROR, "ERROR ! : Error in creating final output file $output_file as $temp_file was not created");
}

###############
# SUBROUTINES #
###############

# Description   : Add snps_per_gene column to the merged table 
# Parameters    : cspg = Column number for snps_per_gene in the merged table
#				  ifile = Path to temporary merged table with all the other information
#				  ofile = Path to final output merged table
#		  		  snps_count = Reference to the hash of count of snps per gene
# Returns       : NA							    
# Modifications :

sub addSnpsPerGene {
	my ($cspg, $ifile, $ofile, $snps_count) = @_;
	my @snps_for_gene = ();
	my @genes = ();
	my @row = ();
	my ($line, $g);
	my $col_gn = $cspg - 5;
	print STDERR "col_gn = $col_gn\n";
	# Loop through each row of the merged table
	open(IN, "< $ifile") or die "Cant open $temp_file\n";
	open(OUT, "> $ofile") or die "Cant open\n";
	while($line = <IN>) {
		chomp($line);
		if($line =~ /^molecule/) {
			print OUT "$line\n";
			next;
		}
		@row = split(/\t/, $line);
		# No need to count snps in intergenic region
		if($row[$col_gn] eq "intergenic") {
			$row[$cspg] = "NA";
			print OUT join("\t", @row);
			print OUT "\n";
			next;
		}
		@snps_for_gene = ();
		# In case of overlapping genes, display the count for both the genes at that SNP position
		@genes = split("/",$row[$col_gn]);
		foreach $g (@genes) {
			if(exists($snps_count->{$row[0]}{$g})) {
				push(@snps_for_gene,$snps_count->{$row[0]}{$g});
			} else {
				printLogMsg($ERROR, "ERROR! : Could not find SNPs per gene for ".$row[0]." and $g");
			}
		}
		$row[$cspg] = join("/", @snps_for_gene);
		print OUT join("\t", @row);
		print OUT "\n";
	}
}

####################################################################################################################################################

# Description   : Parse SNP panel files list to obtain unique SNPs  
# Parameters    : snp_list = path to list of SNP panel files
# Returns       : Reference to a hash of unique SNPs. Key = Molecule name to which the SNP position belongs to, SNP position
#						      Value = Last column in the SNP panel file w.r.t. that SNP position							    
# Modifications :

sub parseSnpPanels {
	my ($snp_file) = @_; 
        my %snp_data = ();
	my @snp_row = ();
	my ($file, $line);
#	my ($temp_k1, $temp_k2);
        open(SL,"< $snp_file") or printLogMsg($ERROR, "ERROR! : Could not open $snp_file file for reading.Reason : $!");
        # Loop through the list of SNP panel files
	while($file = <SL>) {
                chomp($file);
                next if($file =~ /^\s*$/);
                next if($file =~ /^#/);
		# For each SNP panel file, get molecule and SNP position to create a hash of unique SNPs for each molecule
                if(-e $file) {
			printLogMsg($DEBUG, "INFO :: Parsing SNP panel file $file");
                        open(FR,"< $file") or printLogMsg($ERROR, "ERROR! : Could not open $file file for reading.Reason : $!");
                        readline(FR);
                        while(<FR>) {
                                $line = $_; 
                                chomp($line);
                                next if ($line =~ /^\s*$/);
                                next if ($line =~ /^#/);
                                @snp_row = split(/\t/, $line); 
				$snp_row[0] =~ s/^\s+|\s+$//;	
				$snp_row[1] =~ s/^\s+|\s+$//;
				# Store the last column which should be properties or additional info about the SNP position 
				if(exists($snp_data{$snp_row[0]}{$snp_row[1]})) {
					$snp_data{$snp_row[0]}{$snp_row[1]} .= ($#snp_row > 1) ? ("|".$snp_row[$#snp_row]) : ""; 	
				} else {
                                	$snp_data{$snp_row[0]}{$snp_row[1]} = ($#snp_row > 1) ? $snp_row[$#snp_row] : "";
				}
                        }
			close(FR);    
                }
        }
#	foreach $temp_k1 (keys %snp_data) {
#		foreach $temp_k2 (sort {$a <=> $b} keys %{$snp_data{$temp_k1}}) {
#			print "$temp_k1\t\t$temp_k2\t\t$snp_data{$temp_k1}{$temp_k2}\n";
#		}
#	}
	close(SL);
        return (\%snp_data);
}

####################################################################################################################################################

# Description   : Obtain type of SNP - SYN or NSYN, query codon and query amino acid  
# Parameters    : pcodon = relative position in codon of the SNP 
#		  rcodon = reference codon sequence
#		  refnuc = reference base
#		  gstrand = strand on which the gene is located
#		  snptag = SNP molecule and reference position
#		  data = reference to a hash of all SNPs in all the queries
# Returns       : qstr = query codon
#		  qaastr = query amino acid
#		  snpstat = SNP type - SYN or NSYN						      							    
# Modifications :

sub getSnpType {
	my ($pcodon, $rcodon, $refnuc, $gstrand, $snptag, $data) = @_;
	my $qstr = "NA";
	my $qaastr = "NA";
	my $snpstat = "NA";
	my ($queryaa,$stat);
	my (@qcodon,@qaa,@syn);
	# Translate the reference codon to reference amino acid
	my $raa = translateCodon($rcodon);
	# Check if SNP information exists at the SNP position
	if(exists($data->{$snptag})) {
		$qstr = $rcodon;
		$qaastr = $raa;
		$snpstat = "SYN";
		# For all the queries 
		foreach my $q (keys %{$data->{$snptag}}) {
			# For all the bases found at that SNP position for each query
			foreach my $suniq (keys %{$data->{$snptag}{$q}}) {
				next if($suniq =~ /No Hit|indel/); 
				next if($suniq eq $refnuc);
				# Complement the query base if gene is on negative strand
				$suniq =~ tr/ACGT/TGCA/ if ($gstrand == -1);
				# Obtain query codon by replacing query base in the reference codon
				substr($rcodon,($pcodon - 1),1) = $suniq;
				push(@qcodon,$rcodon);
				# Translate query codon to find query amino acid
				$queryaa = translateCodon($rcodon);
				push(@qaa, $queryaa);
				$stat = ($queryaa eq $refaa) ? "SYN": "NSYN";
				push(@syn, $stat);
			}
		}
		if(@qcodon) {
			$qstr = join("/",keys %{{map {$_ => 1}@qcodon}});
			$qaastr = join("/",keys %{{map {$_ => 1}@qaa}});
			$snpstat = join("/",keys %{{map {$_ => 1}@syn}});
		}
	}
	return($qstr, $qaastr, $snpstat);
}

####################################################################################################################################################

# Description   : Translate codon to find coresponding amino acid 
# Parameters    : codon = sequence of the codon
# Returns       : ret = corresponding amino acid or unknown if codon is unusual
# Modifications :

sub translateCodon {
	my ($codon) = @_;
	my $ret = "Unknown";
	# Translation table
	my %aa = (
		'AAA'=>'K',
                'AAC'=>'N',
                'AAG'=>'K',
                'AAT'=>'N',
                'ACA'=>'T',
                'ACC'=>'T',
                'ACG'=>'T',
                'ACT'=>'T',
                'AGA'=>'R',
                'AGC'=>'S',
                'AGG'=>'R',
                'AGT'=>'S',
                'ATA'=>'I',
                'ATC'=>'I',
                'ATG'=>'M',
                'ATT'=>'I',
                'CAA'=>'Q',
                'CAC'=>'H',
                'CAG'=>'Q',
                'CAT'=>'H',
                'CCA'=>'P',
                'CCC'=>'P',
                'CCG'=>'P',
                'CCT'=>'P',
                'CGA'=>'R',
                'CGC'=>'R',
                'CGG'=>'R',
                'CGT'=>'R',
                'CTA'=>'L',
                'CTC'=>'L',
                'CTG'=>'L',
                'CTT'=>'L',
                'GAA'=>'E',
                'GAC'=>'D',
                'GAG'=>'E',
                'GAT'=>'D',
                'GCA'=>'A',
                'GCC'=>'A',
                'GCG'=>'A',
                'GCT'=>'A',
                'GGA'=>'G',
                'GGC'=>'G',
                'GGG'=>'G',
                'GGT'=>'G',
                'GTA'=>'V',
                'GTC'=>'V',
                'GTG'=>'V',
                'GTT'=>'V',
                'TAA'=>'Stop',
                'TAC'=>'Y',
                'TAG'=>'Stop',
                'TAT'=>'Y',
		'TCA'=>'S',
                'TCC'=>'S',
                'TCG'=>'S',
                'TCT'=>'S',
                'TGA'=>'Stop',
                'TGC'=>'C',
                'TGG'=>'W',
                'TGT'=>'C',
                'TTA'=>'L',
                'TTC'=>'F',
                'TTG'=>'L',
                'TTT'=>'F'
		);
	if (($codon =~ /N/) || ($codon =~ /X/) || (length($codon) != 3)){
		return $ret;
	}	
	if (exists($aa{$codon})){
		return $aa{$codon};
	} else {
		return $ret;
	}
}

####################################################################################################################################################

# Description   : Extract a sub-sequence from a hash of sequences given boundary values 
# Parameters    : seq_cache = Reference to a hash of reference genome sequences
#		  parent = molecule name for which sub-sequence is required
#		  start = start of the subsequence
#		  stop = stop for the subsequence
#		  strand = strand for which subsequnce is required
# Returns       : subseq = required subsequence
# Modifications :

sub getSeqFromParent {
	my ($seq_cache, $parent,$start, $stop, $strand) = @_;
	my $seq = "";
	printLogMsg($ERROR, "$parent start $start should be always less than stop $stop") if( $start > $stop );
	# Check if the parent sequence exist 
	if(exists($seq_cache->{$parent})) {
		$seq = $seq_cache->{$parent}; 
	} else {
		printLogMsg($ERROR, "Unable to find FASTA sequence for $parent");
	}
	my $length = ( $stop - $start ) + 1;
	my $subseq = substr( $seq, $start - 1, $length );
	# Reverse complement the subsequence if it is on negative strand
	$subseq = reverseComplement( $subseq ) if( $strand == -1 );
#	print "$subseq from $parent at $start, $stop, $strand\n";
	return ($subseq);
}

####################################################################################################################################################

# Description   : Get reverse complement of a sequence 
# Parameters    : seq = input sequence to be reverse complemented
# Returns       : rev = reverse complemented sequence
# Modifications :

sub reverseComplement {
	my ($seq) = @_;
	my $rev = reverse $seq;
	$rev =~ tr/ATCGRYWSKMVBHDatcgrywskmvbhd/TAGCYRWSMKBVDHtagcyrwsmkbvdh/;
	return($rev);
}

####################################################################################################################################################

# Description   : Create the header line for the output merged table
# Parameters    : queries = Reference to an array of query names 
# Returns       : header = merged table header string
# Modifications :

sub makeHeader {
	my ($queries) = @_;
	my $header = join("\t", ('molecule','refpos','syn?','refbase',@{$queries},'gene_name','product','gene_start','gene_stop','gene_length','snps_per_gene','pos_in_gene','ref_codon','ref_aa','query_codon','query_aa'));
	map { $header .= "\tnum_hits:$_" } @{$queries};
	$header .= "\tproperties";
#	print "$header\n";
	return($header);
}

####################################################################################################################################################

# Description   : Parse raw BLAST results files list to obtain SNP information in query genomes 
# Parameters    : blast_list = Path to a list of BLAST results files
#		  flanking bases = Number of flanking bases used to extract the SNP base from reference genome
#		  qmap = Reference to hash of query genome contig names. Used to map query genome hits back to query name
# Returns       : data = Hash reference containing all SNP information about the hits and queries
# Modifications :

sub parseRawBlast {
	my ($blast_list, $flanking_base, $qmap) = @_;
	my %data = ();
#	my ($temp_k1, $temp_k2, $temp_k3);
	my @match = ();
	my ($file, $query, $in, $result, $snp_pos, $hit, $temp, $hsp, $qstart, $offset, $sseq, $base);
	my $qideal_len = (2 * $flanking_base) + 1;
	open(BL, "< $blast_list") or printLogMsg($ERROR, "ERROR! : Could not open $blast_list file for reading.Reason : $!");
	# Loop through BLAST list file
	while($file = <BL>) {
		chomp($file);
		next if ($file =~ /^\s*$/);
		# Foreach BLAST file
		if (-e $file) {
			printLogMsg($DEBUG, "INFO :: Parsing BLAST file $file");
			$in = "";
			# BLAST file object
			$in = new Bio::SearchIO(-format => 'blast', -file => $file);
			# Foreach SNP position
			while($result = $in->next_result) {
				$snp_pos = $result->query_name;
				# Foreach BLAST hit - there should be one for each query genome
				while($hit = $result->next_hit ) {
					$temp = $hit->name;
					$query = "";
					$temp =~ s/\|/\\\|/g;
					# Query name may be contig names for the genome. So they need to be mapped back to the actual query genome name
					if((exists($qmap->{$hit->name})) || (grep {/$temp/} keys %$qmap)) {
						@match = ();
						$query = $qmap->{$hit->name} if (exists($qmap->{$hit->name}));
						@match = grep { $_ =~ /$temp/;} keys %$qmap if(!length($query));
						$query = $qmap->{$match[0]} if($#match > -1);
					} else {
						printLogMsg($ERROR, "ERROR! : Invalid BLAST hit $temp. Could not find in contigs of the query genomes");
					}
					# For each hsp - multiple hits for each query genome
					while( $hsp = $hit->next_hsp ) {
						$qstart = $hsp->start('query');
						$offset = $flanking_base - ($qstart - 1);
						$sseq = $hsp->hit_string;
						$base = uc(substr( $sseq, $offset, 1 ));
						if($base eq "") {
							$base = "No Hit";
						} elsif ($base eq "-") {
							$base = "indel";
						}
						if(!exists($data{$snp_pos}{$query}{$base})) {
							$data{$snp_pos}{$query}{$base} = 0; 
						} 
						$data{$snp_pos}{$query}{$base}++ ;
					}
				}
			}			
		} else {
			printLogMsg($ERROR, "ERROR! : BLAST file $file does not exist");
		}
	}
#	foreach $temp_k1 (keys %data) {
#			print "$temp_k1\t";
#			my $prev = "";
#		foreach $temp_k2 (keys %{$data{$temp_k1}}) {
#				print "$temp_k2\t";
#			foreach $temp_k3 (keys %{$data{$temp_k1}{$temp_k2}}) {
#				print "$temp_k3" if $temp_k3 ne $prev;
#				$prev = $temp_k3;
#				print "\n";
#			}
#		}
#		print "\n";
#	}
	return(\%data); 
}


####################################################################################################################################################

# Description   : Parse reference GenBank files list to obtain the gene information for each reference genome in the list
# Parameters    : gbk_files = Path to the list of files of reference GenBank genomes
# Returns       : retval = Reference to a hash of gene information. Key = Molecule name, gene name, tags like start, end, length, strand, gene, pseudo, product
#								    Value = Annotation of each gene
#		  seq_cache = Reference to hash of reference genome sequences. Key = Molecule name (display id from GenBank file)
#									       Value = FASTA sequence
# Modifications :

sub parseGenbankFiles {
	my ($gbk_files) = @_;
	my %retval = ();
	my %seq_cache = ();
	my @lts = ();
#	my ($temp_k1, $temp_k2);
	my ($file, $locus_tag, $seq_obj, $accession, $feat, $type, $product);
	open(GL, "< $gbk_files") or printLogMsg($ERROR, "ERROR! : Could not open $gbk_files file for reading.Reason : $!");
	# Loop through each GenBank file
	while($file = <GL>) {
		chomp($file);
		next if ($file =~ /^\s*$/);
		next if ($file =~ /^#/);
		printLogMsg($DEBUG, "INFO :: Parsing GenBank file $file");
		if (-e $file) {
		# Object of GenBank file sequence
			$seq_obj = new Bio::SeqIO( -file => $file )->next_seq();
			$accession = $seq_obj->display_id;
			# Create a hash of each GenBank sequence
			$seq_cache{$accession} = $seq_obj->seq();
			# Foreach sequence feature like gene, CDS
			foreach $feat ($seq_obj->get_SeqFeatures) {
				next if( $feat->primary_tag eq 'source' );				
				@lts = ();
				$type = "";
				$product = "No product";
				if( $feat->primary_tag eq 'gene' ) {
					if( $feat->has_tag('locus_tag') ) {
						@lts = $feat->get_tag_values( 'locus_tag' );
						printLogMsg($ERROR, "ERROR! : Could not parse locus_tag from feature in reference genbank file") unless( @lts );
						printLogMsg($ERROR, "ERROR! : Found more than one locus_tag for feature [@lts] in reference genbank file") if( @lts > 1 );
						$locus_tag = $lts[0];
						printLogMsg($ERROR, "ERROR! : Found multiple genes with same locus_tag: $locus_tag for $accession") if( exists( $retval{$accession}{$locus_tag}{'gene'}));
						$retval{$accession}{$locus_tag}{'start'} = $feat->location->strand == 1 ? $feat->location->start : $feat->location->end;
						$retval{$accession}{$locus_tag}{'end'} = $feat->location->strand == 1 ? $feat->location->end : $feat->location->start;
						$retval{$accession}{$locus_tag}{'length'} = $feat->location->end - $feat->location->start + 1;
						$retval{$accession}{$locus_tag}{'strand'} = $feat->location->strand;
						$retval{$accession}{$locus_tag}{'gene'} = 1;
						$retval{$accession}{$locus_tag}{'pseudo'} = 1 if( $feat->has_tag('pseudo') );
						$retval{$accession}{$locus_tag}{'product'} = $product;
						printLogMsg($ERROR, "ERROR! : Start was greater than stop. [$locus_tag]") if( $feat->location->start > $feat->location->end );
						my $tree_tag = $accession."|".$locus_tag;
						$iTree->addInterval( $tree_tag, $feat->location->start, $feat->location->end );
					} else {
						printLogMsg($ERROR, "ERROR! : Gene feature does not have locus tag");
					}
				} else {
					if( $feat->has_tag('locus_tag') ) {
						@lts = $feat->get_tag_values( 'locus_tag' );
						printLogMsg($ERROR, "ERROR! : Could not parse locus_tag from feature in reference genbank file") unless( @lts );
						printLogMsg($ERROR, "ERROR! : Found more than one locus_tag for feature [@lts] in reference genbank file") if( @lts > 1 );
						$locus_tag = $lts[0];
						$type = $feat->primary_tag;
						printLogMsg($ERROR, "ERROR! : Found multiple $type with same locus_tag: $locus_tag for $accession") if( exists( $retval{$accession}{$locus_tag}{$type} ));
#					if( $feat->location->isa('Bio::Location::SplitLocationI') && $type eq 'CDS' ) {
#						for my $location ( $feat->location->sub_Location ) {
#							$retval{$accession}{$locus_tag}{'start'} = $location->start;
#							$retval{$accession}{$locus_tag}{'end'} = $location->end;
#						}
#					}
						if( $feat->has_tag('product') ) {
							$product = join(", ", $feat->get_tag_values( 'product' ) );
						} elsif($feat->has_tag('note')) {
							$product = join(", ", $feat->get_tag_values( 'note' ) );
						}
						$retval{$accession}{$locus_tag}{$type} = 1; 
					}
					$retval{$accession}{$locus_tag}{'product'} = $product if(exists($retval{$accession}{$locus_tag}));
				}
			}
		} else {
			printLogMsg($ERROR, "ERROR! : GenBank file $file does not exist");
		}
	}
#       foreach $temp_k1 (keys %retval) {
#                       print "$temp_k1\t";
#               foreach $temp_k2 (sort keys %{$retval{$temp_k1}}) {
#                	print "$temp_k2\t $retval{$temp_k1}{$temp_k2}{'start'}\t$retval{$temp_k1}{$temp_k2}{'end'}\t$retval{$temp_k1}{$temp_k2}{'length'}\t$retval{$temp_k1}{$temp_k2}{'strand'}\t$retval{$temp_k1}{$temp_k2}{'product'}";
#                               print "\n";
#               }
#       }	
	return(\%retval, \%seq_cache);
}

####################################################################################################################################################

# Description   : Create a map file to map name of query genomes to the names of the contigs belonging to the query genome  
# Parameters    : query_list = Path to the list 
# Returns       : retval = Reference to a hash of contig names of the query genomes. Key = Contig header
#										     Value = Basename of the query genome file which will be used 
#											     in the final output table as column for the query
#		  queries = Reference to an array of query genome names 
# Modifications :

sub mapQueryList {
	my ($query_list) = @_;
	my %retval = ();
	my @queries = () ;
	my ($file, $file_base, $file_dir, $file_ext);
#	my $temp_k1;
	open(QL, "< $query_list") or printLogMsg($ERROR, "Could not open $query_list file for reading. Reason : $!");
	# Loop through query genome list file
	while($file = <QL>) {
		chomp($file);
		next if ($file =~ /^\s*$/);
		next if ($file =~ /^#/);
		if (-e $file) {
			printLogMsg($DEBUG, "INFO :: Parsing query genome file $file");
			($file_base,$file_dir,$file_ext) = fileparse($file,qr/\.[^.]*/);
			open(QF, "< $file") or printLogMsg($ERROR, "Could not open $file file for reading. Reason : $!");
			while(<QF>) {
				chomp($_);
				if($_ =~ /^>(\S+)/) {
					printLogMsg($ERROR, "Couldn't parse header from defline from fasta file: $_") unless($1);	
					# Hash of contig names and name of the query from the filename
					$retval{$1} = $file_base;
				}		
			}
			close(QF);
			push @queries, $file_base;
		} else {
			printLogMsg($ERROR, "Specified $file file does not exist.");
		}
	}
	close(QL);
#	print "@queries\n";
#	foreach $temp_k1 (sort keys %retval) {
#		print "$retval{$temp_k1}\t\t$temp_k1\n";
#	}
	return(\%retval, \@queries);
}

####################################################################################################################################################

# Description   : Used to check the correctness of the command line arguments passed to the script. The script exits if required arguments are missing. 
# Parameters    : NA
# Returns       : NA
# Modifications :

sub checkCmdLineArgs {
	my @required = qw(raw_blast_list ref_genbank snp_positions query_list output_dir);
	foreach my $option(@required) {
		if(!defined($cmdLineArgs{$option})) {
			printLogMsg($ERROR,"ERROR! : Required option $option not passed");
		}
	}
	if((!exists($cmdLineArgs{'flanking_bases'})) || (!defined($cmdLineArgs{'flanking_bases'}))) {
		$cmdLineArgs{'flanking_bases'} = 20;
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

snp_verify.pl - Script to parse BLAST results and annotate verified SNPs to output merged table. 

=head1 SYNOPSIS

USAGE : perl snp_verify.pl --b <blast_results.list> --r <reference_genbank_file.list> --s <snp_panels.list> --q <query_genome.list> --o <output_directory> [ --f <flanking_bases> --l <log_file> ]

	parameters in [] are optional

=head1 OPTIONS

	--b <blast_results.list> = Path to list file consisting of paths to raw BLAST results. Mandatory

	--r <reference_genbank_file.list> = Path to list file consisting of paths to reference genome GenBank files. Mandatory

	--s <snp_panels.list> = Path to list of files of predicted SNP positions (SNP panel) in the reference genome(s). Mandatory

	--q <query_genome.list> = Path to list file consisting of paths to query genome FASTA files. Mandatory

	--o <output_directory> = /path/to/output_dir. Mandatory

	--f <flanking_bases> = Number of bases on *each side* of the SNP positions in the reference genome. Optional. Default is 20.

	--l <log_file> = /path/to/log_file.log. Optional 

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
