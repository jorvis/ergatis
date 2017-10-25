#!/usr/bin/env perl -w

#########################################################################################
#											#
# Name	      	: create_comparative_summary_report.pl					#
# Version     	: 2.0									#
# Project     	: CloVR Comparative Genomics Pipeline					#
# Description 	: Create summary report of the genomes analysed in the comparative	#
#		  genomics pipeline							#
# Author      	: Sonia Agrawal								#
# Date        	: April 25, 2013							#
# Lasted updated: May 7, 2013								#
#########################################################################################

use strict;
# Extended processing of command line options
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
# Prints usage message from embedded pod documentation
use Pod::Usage;
# Include packages here 
use Bio::SeqIO;
use Math::Round;

###########
# GLOBALS #
###########
my %cmdLineArgs = ();
# Log file handle;
my $logfh;
my ($file,$org, $core_len, $total_cogs, $l, $loc, $cds, $ucds, $len);
my %species = ();
my ($ERROR, $WARN, $DEBUG) = (1, 2, 3);

################
# MAIN PROGRAM #
################
GetOptions(\%cmdLineArgs,
	   'input_gbk_list|i=s',
	   'phylomark_fasta|p=s',
	   'summary_report|r=s',
	   'mugsycog_raw|m=s',
	   'log|l=s', 
	   'help|h'
	  ) or pod2usage();

pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDERR} ) if ($cmdLineArgs{'help'});

checkCmdLineArgs();
open(GL, "< $cmdLineArgs{'input_gbk_list'}") or printLogMsg($ERROR, "ERROR! : Could not open file $cmdLineArgs{'input_gbk_list'} for reading. Reason: $!");
while($file = <GL>) {
	chomp($file);
	next if($file =~ /^\s*$/);
	if (-e $file) {
		parseGenbankFile($file,\%species);
	} else {
		printLogMsg($ERROR, "ERROR! : GenBank file $file does not exist");
	}
}
close(GL);

$core_len = coreGenome($cmdLineArgs{'phylomark_fasta'});

$total_cogs = parseMugsyCog($cmdLineArgs{'mugsycog_raw'}, \%species);

open(SR, "> $cmdLineArgs{'summary_report'}") or printLogMsg($ERROR, "ERROR! : Could not open file $cmdLineArgs{'summary_report'} for writing. Reason: $!");
print SR "************************* SUMMARY REPORT FOR COMPARATIVE GENOMICS PIPELINE FOR ".keys(%species)." INPUT GENOMES *************************\n\n";
print SR "Locus Id\tSpecies\tStrain\tGenome Length(in Mbp)\tNumber of scaffolds/contigs\tNumber of CDS\tNumber of unique CDS\n";
print SR "----------------------------------------------------------------------------------------------------------------------------\n";

foreach $org (keys %species) {
	$cds = 0;
	$ucds = 0;
	$len = 0;
	$l = (keys %{$species{$org}})[-1];
	if(keys(%{$species{$org}}) > 1) {
		print SR "Contigs\t";
	} else {
		print SR "$l\t";
	}
	print SR "$org\t$species{$org}{$l}{'strain'}\t";
	foreach $loc (keys %{$species{$org}}) {
		$len += $species{$org}{$loc}{'length'};
		$cds += $species{$org}{$loc}{'cds'};	
		$ucds += $species{$org}{$loc}{'uniq_cds'};
	}
	print SR nearest(0.01, ($len/1000000))."\t".keys(%{$species{$org}})."\t$cds\t$ucds\n";
}

print SR "\n\n";
print SR "Core genome length (in Kbp) : ".nearest(0.01, ($core_len/1000))."\n\n";
print SR "Number of core clusters : $total_cogs\n\n\n\n";

print SR "References:\n";
print SR "---------\n\n";
print SR "1. CloVR\t\tAngiuoli S.V., Matalka M., Gussman A., Galens K., Vangala M., Riley D., Arze C., White J.R., White O., and Fricke W.F. (2011). CloVR: A virtual machine for automated and portable sequence analysis from the desktop using cloud computing. BMC Bioinformatics, 12(1):356. PMID: 21878105\n\n";
print SR "2. Mugsy\t\tAngiuoli S.V., Salzberg S.L. (2011). Mugsy: fast multiple alignment of closely related whole genomes. Bioinformatics, 27(3):334-42. doi: 10.1093/bioinformatics/btq665. PMID: 21148543\n\n";
print SR "3. Pan-genome annotation\t\tAngiuoli S.V., Dunning Hotopp J.C., Salzberg S.L., Tettelin H. (2011). Improving pan-genome annotation using whole genome multiple alignment. BMC Bioinformatics, 12:272. doi: 10.1186/1471-2105-12-272. PMID: 21718539\n\n";

print SR "4. Phylomark, phylogenetic tree generation\t\tSahl J.W., Matalka M.N., Rasko D.A. (2012). Phylomark, a tool to identify conserved phylogenetic markers from whole-genome alignments. Appl Environ Microbiol, 78(14):4884-92. doi: 10.1128/AEM.00929-12. PMID: 22582056\n\n";

close(SR);

###############
# SUBROUTINES #
###############

# Description   : Used to check the correctness of the command line arguments passed to the script. The script exits if required arguments are missing. 
# Parameters    : NA
# Returns       : NA
# Modifications :

sub parseGenbankFile {
	my ($filename, $spec) = @_;
	my ($feat, $organism, $cds);
	my @strain = ();
#	my $contig = 0;
	my $len;
	my $locus;
	my $flag = 0;
	printLogMsg($DEBUG, "INFO :: Parsing GenBank file $filename");
	# Object of GenBank file sequence
	my $seqio_object = Bio::SeqIO->new(-file => $filename);
	while(my $seq_object = $seqio_object->next_seq()) {
		$cds = 0;
		$locus = $seq_object->display_id;
		$len = length($seq_object->seq());
		$organism = $seq_object->species->node_name;
		foreach $feat ($seq_object->get_SeqFeatures) {
			if($feat->primary_tag eq 'source') {
				if($feat->has_tag("strain")) {
					push @strain, $feat->get_tag_values("strain");
					$flag = 1;
				} else {
					$strain[0] = "-" if($flag == 0);
				}
			}
			if($feat->primary_tag eq 'CDS') {
				$cds++;
			}
		}
		$spec->{$organism}{$locus}{'length'} = $len;
		$spec->{$organism}{$locus}{'cds'} = $cds;
		$spec->{$organism}{$locus}{'strain'} = $strain[0];	
		$spec->{$organism}{$locus}{'uniq_cds'} = 0;
	}
}

####################################################################################################################################################

# Description   : Used to check the correctness of the command line arguments passed to the script. The script exits if required arguments are missing. 
# Parameters    : NA
# Returns       : NA
# Modifications :

sub coreGenome {
        my ($filename) = @_; 
        my ($head, $seq, $line, $prev, $cores);
        my %genome = (); 
        open(FR, "< $filename") or printLogMsg($ERROR, "Could not open file $filename for reading. Reason: $!");
        while($line = <FR>) {
                chomp($line);
                next if ($line =~ /^\s*$/);
                if($line =~ /^>(.+)/) {
                        if($seq) {
                                $genome{$head} = length($seq);  
                        }
                        $head = $1; 
                        $head =~ s/\s+$//;    
                        $seq = "";    
                } else {
                        $seq .= $line;
                }
        }
        $genome{$head} = length($seq);
        close(FR);
	foreach $cores (keys %genome) {
		if(!$prev) {
			$prev = $genome{$cores};
			next;
		}
		if($genome{$cores} != $prev) {
			printLogMsg($ERROR, "Core genome length for $cores does not match with other genomes. File $filename doesnot contain correct information about core genome");
		}
		$prev = $genome{$cores};
	}
        return($prev);
}

####################################################################################################################################################

# Description   : Used to check the correctness of the command line arguments passed to the script. The script exits if required arguments are missing. 
# Parameters    : NA
# Returns       : NA
# Modifications :

sub parseMugsyCog {
	my ($filename, $spec) = @_;
	my ($line, $num_genomes, $s, $str, $org);
	my $core_cogs = 0;
	my @temp = ();
	open(FR, "< $filename") or printLogMsg($ERROR, "Could not open file $filename for reading. Reason: $!");
	$num_genomes = keys(%$spec);
	while($line = <FR>) {
		chomp($line);
		if($line =~ /^>/) {
			if($line =~ /num_seqs=(\d+)/) {
				$core_cogs++ if($1 == $num_genomes);
			}
		} elsif($line =~ /^#SINGLETON\s(\S+)/) {
			($org, $str) = split(/\|\|\|/, $1, 2);
			foreach $s (keys %$spec) {
				if(exists($spec->{$s}{$org})) {
					if(exists($spec->{$s}{$org}{'uniq_cds'})) {
						$spec->{$s}{$org}{'uniq_cds'}++;
					} else {
						$spec->{$s}{$org}{'uniq_cds'} = 1;
					}
				}  
			}
		}
	}
	close(FR);
	return($core_cogs);
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
	my @required = qw(input_gbk_list phylomark_fasta summary_report);
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
