#! /usr/local/bin/perl -w


use strict;
use lib "/home/jravel/lib/";
use BeginPerlBioinfo;
use Getopt::Long;
$| = 1;

my $coords = "/home/jravel/src/SNP/GBA_com.coords";
my $genome = "/home/jravel/src/SNP/GBA6615.2con";
my $tag;
my $filename;
my $list = "SNP_location.txt";
my %orf;
my $asmbl_id;
my $posR;
my $posA;
my $end3;
my $end5;
my $end5_3;
my $orf_name;
my $lenA;
my $dist;
my $mod;
my $frame;
my $aa_seq;
my @asmbl;
my $aa_seq2;
my $aseq;

my $help;
my $HELPTEXT = qq~


~;
    


##############################################################################
### OPTION TO RUN THE PROGRAM
##############################################################################

Getopt::Long::config("no_ignore_case");

my $result = GetOptions ("h|help"     => \$help,
			 "F=s"        => \$filename,
			 "T=s"        => \$tag);
	

if ($result == 0) {
    die ("Command line parsing failed\n")
    }

if ($help) {

    print $HELPTEXT, "\n";
    exit;
}

if (! defined $filename) {
    print STDERR "You must enter an input file (-F)\n";
    exit(0);
}

if (! defined $tag) {
    print STDERR "You must enter a folder tag (-T)\n";
    exit(0);
}

my $directory = "../" . $tag ."_asmbl_seq.dir";

#############################################################################
### CREATE LOOK UP FOR COORDINATES
#############################################################################


my @coords_data = get_file_data($coords);

foreach my $line (@coords_data) {

    my @tp = split("::", $line);

    
    $orf{$tp[0]} = "$tp[1]::$tp[2]";

    
}

#############################################################################
###  EXTRACT SEQUENCE FROM 1con
#############################################################################

my @genome = get_file_data($genome);

my $gseq = extract_sequence_from_fasta_data(@genome);


#############################################################################
### READ INPUT FILE : SNP LIST FROM SNP_locaton.txt
###
### SNP_id  REF LENGTH      REF_POS SNP     QUE_POS ORF_ID  COM_NAME
### >ID478revcom-1  (26240bp)       5211344 t-->c   725     ORF01755        DHH subfamily 1 protein
### >ID478revcom-2  (26240bp)       5212286 t-->g   1667    ORF01755        DHH subfamily 1 protein
### >ID478revcom-3  (26240bp)       5213043 c-->t   2424    ORF01754        conserved hypothetical protein
### >ID478revcom-4  (26240bp)       5213582 a-->g   2963    Intergenic
### >ID478revcom-5  (26240bp)       5217760 a-->g   7141    ORF01745        conserved hypothetical protein
### >ID478revcom-6  (26240bp)       5218510 a-->g   7891    ORF01744        stage 0 sporulation protein J
### >ID478revcom-7  (26240bp)       5225411 t-->c   14792   ORF01738        jag protein
###
#############################################################################


my @snp = get_file_data($list);

my $i = 0;
my $int = 0;
my $in = 0;

foreach my $line (@snp) {


    print "THIS IS LINE $line\n";

    if ($i == 0) {

	$i++;

	next;
    }

    my @tmp = split(" ", $line);

    print "THIS IS TEMP: $tmp[0]\n$tmp[1]\n$tmp[2]\n";

    
    if ($tmp[5] eq "Intergenic") {
	
	$int++;
	next;
    }elsif ($tmp[5] =~ /ORF/) {


	($asmbl_id) = ($tmp[0] =~ /^\>(ID\w+)-.*/);

	print "THIS IS ASMBL_ID: $asmbl_id \n";


	($lenA) = ($tmp[1] =~ /^\((\d+)bp.*/);

	$posR = $tmp[2];
	print "THIS IS POSR: $posR \n";
	$posA = $tmp[4];
	print "THIS IS POSA: $posA \n";
	$orf_name = $tmp[5];
	my $trans = $tmp[3];

	print "THIS IS SNP: $trans\n";


	$end5_3 = $orf{$orf_name};
	
	my @tp1 = split("::", $end5_3);
	
	$end5 = $tp1[0];
	$end3 = $tp1[1];

	if ($end3 > $end5) {

	    $dist = $posR - $end5 + 1;
	    print "THIS IS DISTANCE: $dist \n";

	    $mod = $dist%3;

	    print "THIS IS MODULO: $mod \n";


	    if ($mod == 0) {

		$aa_seq = quickcut($gseq, $posR-2, $posR);
		
		print "THIS IS AA_SEQ: $aa_seq\n";

		print "THIS IS AA: " . codon2aa($aa_seq) ."\n\n";


		@asmbl = get_file_data("$directory/$asmbl_id");
		$aseq = extract_sequence_from_fasta_data(@asmbl);

		$aa_seq2 = quickcut($aseq, $posA-2, $posA);

		print "THIS IS AA_SEQ_A: $aa_seq2\n";

		print "THIS IS AA_A: " . codon2aa($aa_seq2) . "\n\n";


	    }

	     if ($mod == 2) {

		$aa_seq = quickcut($gseq, $posR-1, $posR+1);
		
		print "THIS IS AA_SEQ: $aa_seq\n";

		print "THIS IS AA: " . codon2aa($aa_seq) ."\n\n";


		@asmbl = get_file_data("$directory/$asmbl_id");
		$aseq = extract_sequence_from_fasta_data(@asmbl);

		$aa_seq2 = quickcut($aseq, $posA-1, $posA+1);

		print "THIS IS AA_SEQ_A: $aa_seq2\n";

		print "THIS IS AA_A: " . codon2aa($aa_seq2) . "\n\n";


	    }

	    if ($mod == 1) {

		$aa_seq = quickcut($gseq, $posR, $posR+2);
		
		print "THIS IS AA_SEQ: $aa_seq\n";

		print "THIS IS AA: " . codon2aa($aa_seq) ."\n\n";


		@asmbl = get_file_data("$directory/$asmbl_id");
		$aseq = extract_sequence_from_fasta_data(@asmbl);

		$aa_seq2 = quickcut($aseq, $posA, $posA+2);

		print "THIS IS AA_SEQ_A: $aa_seq2\n";

		print "THIS IS AA_A: " . codon2aa($aa_seq2) . "\n\n";


	    }

	}elsif ($end3 < $end5) {

	    $dist = $end5 - $posR+ 1;
	    print "THIS IS DISTANCE: $dist \n";

	    $mod = $dist%3;

	    print "THIS IS MODULO: $mod \n";


	    if ($mod == 0) {

		$aa_seq = quickcut($gseq, $posR, $posR+2);
		
		$aa_seq = revcom($aa_seq);
		
		print "THIS IS AA_SEQ: $aa_seq\n";

		print "THIS IS AA: " . codon2aa($aa_seq) ."\n\n";


		@asmbl = get_file_data("$directory/$asmbl_id");
		$aseq = extract_sequence_from_fasta_data(@asmbl);

		$aa_seq2 = quickcut($aseq, $posA, $posA+2);

		$aa_seq2 = revcom($aa_seq2);

		print "THIS IS AA_SEQ_A: $aa_seq2\n";

		print "THIS IS AA_A: " . codon2aa($aa_seq2) . "\n\n";


	    }

	     if ($mod == 2) {

		$aa_seq = quickcut($gseq, $posR-1, $posR+1);
		
		$aa_seq = revcom($aa_seq);

		print "THIS IS AA_SEQ: $aa_seq\n";

		print "THIS IS AA: " . codon2aa($aa_seq) ."\n\n";


		@asmbl = get_file_data("$directory/$asmbl_id");
		$aseq = extract_sequence_from_fasta_data(@asmbl);

		$aa_seq2 = quickcut($aseq, $posA-1, $posA+1);

		$aa_seq2 = revcom($aa_seq2);

		print "THIS IS AA_SEQ_A: $aa_seq2\n";

		print "THIS IS AA_A: " . codon2aa($aa_seq2) . "\n\n";


	    }

	    if ($mod == 1) {

		$aa_seq = quickcut($gseq, $posR-2, $posR);
		
		$aa_seq = revcom($aa_seq);

		print "THIS IS AA_SEQ: $aa_seq\n";

		print "THIS IS AA: " . codon2aa($aa_seq) ."\n\n";


		@asmbl = get_file_data("$directory/$asmbl_id");
		$aseq = extract_sequence_from_fasta_data(@asmbl);

		$aa_seq2 = quickcut($aseq, $posA-2, $posA);
		
		$aa_seq2 = revcom($aa_seq2);

		print "THIS IS AA_SEQ_A: $aa_seq2\n";

		print "THIS IS AA_A: " . codon2aa($aa_seq2) . "\n\n";


	    }
	    
	}
	
    }
}

exit;
