#! /usr/local/bin/perl -w


use strict;
use lib "/home/jravel/lib/";
use BeginPerlBioinfo;
use Getopt::Long;
$| = 1;

my $coords = "/home/jravel/src/SNP/GBX_pXO1_com.coords";
my $genome = "/home/jravel/src/SNP/GBX_01.1con";
my $tag;
my $filename;
my $list = "SNP_final.txt";
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
my @snp;
my @tmp2;
my @tmp;
my $aa_seq2;
my $aseq;
my $record;
my $fh;
my $newline;
my $syn;
my $ctsyn;
my $ctnsyn;
my $aa_r;
my $aa_a;
my $orf_len;
my $codpos;
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

my $directory = $tag ."_asmbl_seq.dir";

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

open($fh, $list);


my $i = 0;
my $int = 0;
my $in = 0;


while ($record = get_next_record4($fh)) {


    @tmp2 = split("\n", $record);

    print "THIS IS TMP2[4]: $tmp2[4]\n";
   
    if ($tmp2[4] =~ /^Intergenic/) {
	
	chomp $record;
	push(@snp, "$record");
	$int++;
	next;

    }elsif ($tmp2[4] =~ /ORF/) {


	@tmp = split(" ", $tmp2[0]);

	($orf_name) = ($tmp2[4] =~ /^(ORF.\d+)\s+.*/);

	print "THIS IS ORF_NAME: $orf_name\n";

	($asmbl_id) = ($tmp[0] =~ /^\>(ID\d+\w+)-.*/);

	print "THIS IS ASMBL_ID: $asmbl_id \n";


	($lenA) = ($tmp[1] =~ /^\((\d+)bp.*/);

	$posR = $tmp[2];
	print "THIS IS POSR: $posR \n";
	$posA = $tmp[4];
	print "THIS IS POSA: $posA \n";
	
	my $trans = $tmp[3];

	print "THIS IS SNP: $trans\n";


	$end5_3 = $orf{$orf_name};
	
	my @tp1 = split("::", $end5_3);
	
	$end5 = $tp1[0];
	$end3 = $tp1[1];

	if ($end3 > $end5) {

	    $orf_len = $end3 - $end5;

	    $dist = $posR - $end5 + 1;
	    print "THIS IS DISTANCE: $dist \n";

	    $mod = $dist%3;

	    print "THIS IS MODULO: $mod \n";


	    if ($mod == 0) {

		$codpos = 3;

		$aa_seq = quickcut($gseq, $posR-2, $posR);
		
		print "THIS IS AA_SEQ: $aa_seq\n";

		$aa_r = codon2aa($aa_seq);

		print "THIS IS AA: " . codon2aa($aa_seq) ."\n\n";


		@asmbl = get_file_data("$directory/$asmbl_id");
		$aseq = extract_sequence_from_fasta_data(@asmbl);

		$aa_seq2 = quickcut($aseq, $posA-2, $posA);

		print "THIS IS AA_SEQ_A: $aa_seq2\n";

		$aa_a = codon2aa($aa_seq2);

		print "THIS IS AA_A: " . codon2aa($aa_seq2) . "\n\n";


	    } elsif ($mod == 2) {

		$codpos = 2;

		$aa_seq = quickcut($gseq, $posR-1, $posR+1);
		
		print "THIS IS AA_SEQ: $aa_seq\n";

		$aa_r = codon2aa($aa_seq);

		print "THIS IS AA: " . codon2aa($aa_seq) ."\n\n";


		@asmbl = get_file_data("$directory/$asmbl_id");
		$aseq = extract_sequence_from_fasta_data(@asmbl);

		$aa_seq2 = quickcut($aseq, $posA-1, $posA+1);

		print "THIS IS AA_SEQ_A: $aa_seq2\n";

		$aa_a = codon2aa($aa_seq2);

		print "THIS IS AA_A: " . codon2aa($aa_seq2) . "\n\n";


	    }elsif ($mod == 1) {
		
		$codpos = 1;

		$aa_seq = quickcut($gseq, $posR, $posR+2);
		
		print "THIS IS AA_SEQ: $aa_seq\n";

		$aa_r = codon2aa($aa_seq);

		print "THIS IS AA: " . codon2aa($aa_seq) ."\n\n";


		@asmbl = get_file_data("$directory/$asmbl_id");
		$aseq = extract_sequence_from_fasta_data(@asmbl);

		$aa_seq2 = quickcut($aseq, $posA, $posA+2);

		print "THIS IS AA_SEQ_A: $aa_seq2\n";

		$aa_a = codon2aa($aa_seq2);

		print "THIS IS AA_A: " . codon2aa($aa_seq2) . "\n\n";


	    }

	}elsif ($end3 < $end5) {

	    $orf_len = $end5 - $end3;

	    $dist = $end5 - $posR+ 1;
	    print "THIS IS DISTANCE: $dist \n";

	    $mod = $dist%3;

	    print "THIS IS MODULO: $mod \n";


	    if ($mod == 0) {

		$codpos = 3;

		$aa_seq = quickcut($gseq, $posR, $posR+2);
		
		$aa_seq = revcom($aa_seq);
		
		print "THIS IS AA_SEQ: $aa_seq\n";
		
		$aa_r = codon2aa($aa_seq);
		
		print "THIS IS AA: " . codon2aa($aa_seq) ."\n\n";


		@asmbl = get_file_data("$directory/$asmbl_id");
		$aseq = extract_sequence_from_fasta_data(@asmbl);

		$aa_seq2 = quickcut($aseq, $posA, $posA+2);

		$aa_seq2 = revcom($aa_seq2);

		print "THIS IS AA_SEQ_A: $aa_seq2\n";

		$aa_a = codon2aa($aa_seq2);

		print "THIS IS AA_A: " . codon2aa($aa_seq2) . "\n\n";


	    }elsif ($mod == 2) {

		$codpos = 2;

		$aa_seq = quickcut($gseq, $posR-1, $posR+1);
		
		$aa_seq = revcom($aa_seq);

		print "THIS IS AA_SEQ: $aa_seq\n";

		$aa_r = codon2aa($aa_seq);

		print "THIS IS AA: " . codon2aa($aa_seq) ."\n\n";


		@asmbl = get_file_data("$directory/$asmbl_id");
		$aseq = extract_sequence_from_fasta_data(@asmbl);

		$aa_seq2 = quickcut($aseq, $posA-1, $posA+1);

		$aa_seq2 = revcom($aa_seq2);

		print "THIS IS AA_SEQ_A: $aa_seq2\n";

		$aa_a = codon2aa($aa_seq2);

		print "THIS IS AA_A: " . codon2aa($aa_seq2) . "\n\n";


	    }elsif ($mod == 1) {

		$codpos = 1;

		$aa_seq = quickcut($gseq, $posR-2, $posR);
		
		$aa_seq = revcom($aa_seq);

		print "THIS IS AA_SEQ: $aa_seq\n";

		$aa_r = codon2aa($aa_seq);

		print "THIS IS AA: $aa_r \n\n";


		@asmbl = get_file_data("$directory/$asmbl_id");
		$aseq = extract_sequence_from_fasta_data(@asmbl);

		$aa_seq2 = quickcut($aseq, $posA-2, $posA);
		
		$aa_seq2 = revcom($aa_seq2);

		print "THIS IS AA_SEQ_A: $aa_seq2\n";
		
		$aa_a = codon2aa($aa_seq2);

		print "THIS IS AA_A: $aa_a \n\n";


	    }
	    
	}

    }
    
	
    if ($aa_r eq $aa_a) {
	    
	$syn = "SYNONYMOUS";
	print $syn , "\n";
	
	$ctsyn++;
    }else  {
	
	$syn = "NON-SYNONYMOUS";
	print $syn , "\n";
	$ctnsyn++;
    }
    
    $dist -= 1;

    $newline = "$syn ($codpos): $aa_seq - $aa_r ---> $aa_seq2 - $aa_a   (from start: $dist bp/$orf_len bp)";

    print "THIS IS NEWLINE: $newline\n\n";
	
    my $zt = 0;

    foreach my $line (@tmp2) {

	if ($zt == 1) {
	    push(@snp, $newline);
	    push(@snp, $line);
	    $zt++;
	    next;
	}else {
	    push(@snp, $line);
	    $zt++;
	    next;
	}
    }
    
}



open(OUT, ">output.final");

foreach my $line (@snp) {
    
    print OUT $line, "\n";;
}


close (OUT);

exit;


sub get_next_record4 {

    my ($fh) = @_;

    my $offset;
    my $record = '';
    my ($save_input_separator) = $/;

    $/ = "//\n";

    $record = <$fh>;
    $/ = $save_input_separator;

    return $record;

}

