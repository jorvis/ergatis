#! /usr/local/bin/perl -w

#######################################################################################################
### A script to cut out dna sequences given two position surrounding a SNP       
### using getTotalCoverage generated files.
###                             
### Input = output from read_from_align.pl (T. Read script)
### Output = multifasta file of sequence surrounding the SNP
###  ID# (size of assembly) position on REF SNP(bp-->bp) position on assembly
###  position on assembly of the sequence (500 bp around SNP).
###  user define outputname with option -O filename
### output = a list of all the headers from the multifasta
###  SNP_Header.txt
### output = a list of alignment 20 bp aournd the SNP with the
### reference genome.
###  SNP_align.txt
###  SNP_quality.txt
###  
### USAGE: fasta_SNP.pl -G genome.fasta -F asemmbly.fasta -I mummer.align -C assembly.contig -T tag 
### -D asmbl_id -V contig.tcov -O outputfilename
###
### Last update 10/30/02 by jravel@tigr.org
###
### 
### WORK IN ROGRESS TO INCORPORATE THE ANALYSIS OF INDELS Started 10/30/02
###
### For more information please read the README.doc file. This file is part of a package.
#######################################################################################################


use strict;
use lib "/home/jravel/lib/";
use sub;
use Getopt::Long;
$| = 1;

##### DECLARE VARIABLE

my ($readfile);  ## MUMMER .align file
my ($outputfile);
my ($length) = "100";
my ($pos1);
my ($pos2);
my ($help);
my ($fasta); ## ASSEMBLY FILE IN FASTA FORMAT
my (@newarray);
my ($entry2);
my ($entry3);
my ($entry4);
my (@header);
my ($header) = "SNP_Header.txt";
my ($genome);  ### REFERENCE GENOME CHANGE HERE TO ENTER AT PROMPT
my ($apos1);
my ($apos2);
my ($gpos1);
my ($gpos2);
my (@align);
my (@galign);
my ($subgseq);
my ($subaseq);
my ($alignfile) = "SNP_align.txt";
my ($contigfile); ### .contig FILE FROM PULLCONTIG WILL NEED TO BE PIPED
my ($qualityfile) = "SNP_quality.txt";
my ($asmbl_id); ### THE ASSEMBLY ID NUMBER FOR CutAsm
my (@quality);
my ($posM);
my (@temp);
my $qual2;
my $sum;
my $comment;
my $len;
my $coverage;
my $val;
my $lval;
my $lval2;
my $pval;
my $elem;
my $val2;
my $star;
my $tcov;
my @coverage;
my $newcov;
my @cov3;
my $tag;
my @indelR;
my @indelQ;
my $indelR = "refins.txt";
my $indelQ = "queryins.txt";
my @snpseq;
my $snpout = "snpqual.txt";

Getopt::Long::config("no_ignore_case");

my $result = GetOptions ("h|help"     => \$help,
			 "G=s"        => \$genome,
			 "F=s"        => \$fasta,
			 "I=s"        => \$readfile,
			 "C=s"        => \$contigfile,
			 "T=s"        => \$tag,
			 "D=s"        => \$asmbl_id,
			 "V=s"        => \$tcov,
			 "O=s"        => \$outputfile);

if ($result == 0) {
    die ("Command line parsing failed\n")
    }
if (! defined $fasta) {
    print STDERR "You must enter a valid fasta file\n";
    exit(0);
}
if (! defined $genome) {
    print STDERR "You must enter a valid genome file\n";
    exit(0);
}
if (! defined $readfile) {
    print STDERR "You must enter a valid input file\n";
    exit(0);
}
if (! defined $contigfile) {
    print STDERR "You must enter a valid contig file\n";
    exit(0);
}

if (! defined $tag) {
    print STDERR "You must specify a tag folder name\n";
    exit(0);
}

if (! defined $asmbl_id) {
    print STDERR "You must enter a valid assembly id\n";
    exit(0);
}
if (! defined $tcov) {
    print STDERR "You must specify a coverage file name\n";
    exit(0);
}
if (! defined $outputfile) {
    print STDERR "You must specify a output file name\n";
    exit(0);
}


#### Get the data from the seqfile and put it into an array ####

## FOR THE ASSEMBLY

my @dna = get_file_data($fasta);

## FOR THE REFERENCE GENOME

my @genome = get_file_data($genome);



#### Extract the sequence data from the array and remove blank space ####

my $dna = extract_sequence_from_fasta_data(@dna);

my $gseq = extract_sequence_from_fasta_data(@genome);

my $length2 = length($dna); ### LENGTH OF THE CONTIG
 
print $length2, "\n";


my $tagname = $tag . "_asmbl_seq.dir"; #### NAME OF THE FOLDER WHERE THE FASTA FILES ARE LOCATED

$fasta =~ s/^$tagname\///;        ### REMOVE THE NAME OF THE FOLDER FROM THE NAME FILE


#### OPEN THE COVERAGE FILE FOR THIS CONTIG AND STORE IT INTO AN ARRAY
#### THIS ARRAY WILL BE LOOKED UP WHEN THE COVERAGE/QUALITY VALUES WILL BE FORMATTED

open(QUAL, $tcov);
	
@coverage = <QUAL>;
close (QUAL);

my $t;
my @at;

for ($t = 0; $t < @coverage; $t++) {

    if ($coverage[$t] =~ /^(\d+)/ && !$at[$1]) {
	$at[$1] = $t;
    }
}

#print $at[91361], "\n";

#print $coverage[$at[91361]], "\n";

#print  @coverage[$at[91356]..$at[91366]], "\n";

#print @covtemp, "\n";



## CALL READ_MUMMER_ALIGN AND PROCESS THE MUMMER.align FILE

#  open(READ, "~/src/SNP/read_mumer_align2.perl $readfile |");   

open (READ, "~/src/SNP/change_align.pl $readfile |");

@temp = <READ>;
close (READ);

print @temp;

#    my @temp = get_file_data($readfile);  ## TO BE DELETED ONCE ABOVE IMPLEMENTED
### OTHERWISE USE read_mummer_align.perl which output to a filexc


    my (@read);
    my $i=1;
    
    foreach my $line (@temp) {
       
	chomp $line;

	@read = split(" ", $line);

	### DOES NOT PROCESS GAPS (INSERTION/DELETION)
	### IF . IN SECOND COLUMN (REF) SKIP
       
	if ($read[1] =~ /-/){
	    $line .= "\t$fasta";
	    push(@indelQ, $line);
	    next;
	}

        ### DOES NOT PROCESS GAPS (INSERTION/DELETION)
	### IF . IN SECOND COLUMN (REF) SKIP

	if ($read[2] =~ /-/){
	    $line .= "\t$fasta";
	    push(@indelR, $line);
	    next;
	}

	### DOES NOT PROCESS AMBIGUOUS BP

	if ($read[2] =~/[mrwsykn]/) {
	    next;
	}

	### LAST LINE SKIP

	if ($read[0] =~ /Total/){
	    next;

	### IF A SNP PROCESS:
    
	}else{
	    
	    ### DETERMINE THE POSITION OF THE SNP COMPARE TO THE ALL ASMBL
	    ### TEST IF ASSEMBLY LESS THAN 500 and CALCULATE THE CUT-OUT ACCORDINGLY
	    
	    if ($read[3] > 500 and $read[3] < ($length2 - 500)){
		$pos1 = $read[3] - 500; ##FOR FASTA FILE 500 bp around the SNP
		$pos2 = $read[3] + 500;
		$comment = "SNP: 501/1001";
	    }elsif ($read[3] < 500 and $read[3] < ($length2 - 500)) {   #### IF THE SNP IS A LESS THAN 500 bp 
		$pos1 = 1;                                              #### FROM THE START OF ASMBL
		$pos2 = $read[3] + 500;
		$comment = "SNP: $read[3]/$pos2";
	    }elsif ($read[3] > 500 and $read[3] > ($length2 - 500)) {   #### IF SNP AT LESS THAN 500 bp
		$pos1 = $read[3] - 500;                                 #### FROM THE END OF ASMBL
		$pos2 = $length2;
		$len = ($length2 - $read[3] + 501);
		$comment = "SNP: 501/$len";
	    }elsif ($read[3] < 500 and $read[3] > ($length2 - 500)) {   ### IF SNP AT LESS THAN 500 bp 
		$pos1 = 1;                                              ### FROM START OR END
		$pos2 = $length2;
		$len = $length2;;
		$comment = "SNP: $read[3]/$len";
	    }elsif ($read[3] < 500 and $length2 < 501) {
		$pos1 = 1;
		$pos2 = $length2;
		$comment = "SNP: $read[3]/$pos2";
	    }
		
	    $apos1 = $read[3] - 20; ## FOR ALIGN FILE 20 bp around the SNP ON ASSEMBLY
	    $apos2 = $read[3] + 20;
	    $gpos1 = $read[0] - 20; ## FOR ALIGN FILE 20 bp around the SNP ON GENOME
	    $gpos2 = $read[0] + 20;

	if ($fasta =~ /revcom/){  ## IF THE FASTA FILE CONTAIN revcom THEN COMPUTE COORDINATE
	    $posM = ($length2 - $read[3] + 1);  ### LENGTH OF ASSEMBLY - SNP pos + 1 (start at 0)

	}else{  ## IF ON FORWARD STRAND JUST KEEP THE POSITION OF THE SNP
	    $posM = $read[3];
        }

##################################################################
####     FORMAT THE OUTPUT
##################################################################
	    

	   
	    print $fasta;
	    print $posM, "\n\n";


	### FORMAT THE FILE TO GET QUALITY VALUES

	    push(@quality, "\n>" . $fasta . "-" . $i . " " .  $read[0] ." ". $read[1] ."-->". $read[2] ." ". $read[3]);
	    push(@snpseq,  ">" . $fasta . "-" . $i . " " .  $read[0] ." ". $read[1] ."-->". $read[2] ." ". $read[3]);

	    ### LOOK UP THE COVERAGE ARRAY AND PULL THE LINE FOR +/- 5 around the SNP

	    @cov3 = @coverage[$at[$posM-5]..$at[$posM+5]];
	    
	    foreach my $covline (@cov3) { ### GO THROUGH EACH LINE
		
		my @cov2 = split(" ", $covline);

		if (!defined $cov2[3]) { $cov2[3] = "NC";}
		
		if (!defined $cov2[4]) { $cov2[4] = "NC";}

		
		if ($cov2[0] == $posM) {  ### TAG THE LINE FOR THE SNP WITH *
		    
		    $pval = $cov2[2];
		    $lval = length($cov2[3]);
		    
		   
		    
		    $newcov = "*" . "$cov2[0] \t $cov2[1] \t $cov2[2] \t $cov2[3]   \t $cov2[4]";
		    push (@quality, $newcov);

		    push (@snpseq, $newcov);

		}else {  ### BASICALLY DO NOTHING HERE BUT ADD THE TAB

		    $newcov = " " . "$cov2[0] \t $cov2[1] \t $cov2[2] \t $cov2[3]   \t $cov2[4]";
		    push (@quality, $newcov);

		}

	    }  
		
	     
	   
	## HEADER

	my $alignhead = ">" . $fasta . "-" . $i." (".$length2."bp) " .  $read[0] ." ". $read[1] ."-->". $read[2] ." ". $read[3]." --- " . $apos1. " to " . $apos2 ." COV: " . $lval . " CB_QVal: ". $pval;;

	push(@galign, $alignhead);

	## GET THE SEQUENCE AROUND THE SNP 20 bp AROUND/ CALL TO QUICKCUT AS A SUBROUTINE IN SUB

	$subgseq = quickcut ($gseq, $gpos1, $gpos2);

	## FORMAT THE ALIGNMENT

	push(@galign, "REF " .  ": " . $subgseq . " " . $gpos2);
	push(@galign, "      |||||||||||||||||||| ||||||||||||||||||||");
	$subaseq = quickcut ($dna, $apos1, $apos2);
	push(@galign, "ASM " .  ": " . $subaseq . " " . $apos2. "\n");
       
	## FORMAT THE FASTA FILE OF SEQUENCE 500 bp AROUND THE SNP

	    ## FORMAT THE HEADER 

	my $entry = ">" . $fasta . "-" . $i." (".$length2."bp) " .  $read[0] ." ". $read[1] ."-->". $read[2] ." ". $read[3]." --- " . $pos1. " to " . $pos2 . " (" . $comment . ") COV: " . $lval . " CB_QVal: ". $pval;
	
	push(@newarray, $entry);

	push(@header, $entry);

            ## GET THE SEQUENCE AROUND THE SNP 500 bp on the ASSEMBLY

	$entry2 = quickcut ($dna, $pos1, $pos2);
	
	   ## FORMAT THE SEQUENCE IN LINES OF 100 bp (CAN BE CHANGED AT $length

	for ( my $pos = 0 ; $pos < length($entry2) ; $pos += $length ) {
        push(@newarray, substr($entry2, $pos, $length));
        }
    

	

	$i++;   ### COUNTER TO GIVE A NUMBER TO THE SNP
    }
    }
	
### PRINT TO FILES

    ### PRINT TO SNP_quality.txt

unless (open(QUALITY, ">>$qualityfile")) {
    print STDERR "Can't open file!!\n\n";
    exit;
}

foreach my $sub (@quality) { 
    print QUALITY $sub, "\n";
}

close (QUALITY);

    ### PRINT TO SNP_align.txt

unless (open(ALIGN, ">>$alignfile")) {
    print STDERR "Can't open file!!\n\n";
    exit;
}

foreach my $sub (@galign) { 
    print ALIGN $sub, "\n";
}

close (ALIGN);

    ### PRINT TO SNP_Header.txt

unless (open(HEAD, ">>$header")) {
    print STDERR "Can't open file!!\n\n";
    exit;
}

foreach my $head (@header) { 
    print HEAD $head, "\n";
}

close (HEAD);

    ### PRINT TO SNP_fasta.txt

unless (open(FILEOUT, ">>$outputfile")) {
    print STDERR "Can't open file!!\n\n";
    exit;
}

foreach my $line2 (@newarray) { 
    print FILEOUT $line2, "\n";
}

close (FILEOUT);

unless (open(SNPOUT, ">>$snpout")) {
    print STDERR "Can't open file!!\n\n";
    exit;
}

foreach my $line5 (@snpseq) { 
    print SNPOUT $line5, "\n";
}

close (FILEOUT);

unless (open(INDELR, ">>$indelR")) {
    print STDERR "Can't open file!!\n\n";
    exit;
}

foreach my $indel (@indelR) { 
    print INDELR $indel, "\n";
}

close (INDELR);

unless (open(INDELQ, ">>$indelQ")) {
    print STDERR "Can't open file!!\n\n";
    exit;
}

foreach my $delin (@indelQ) { 
    print INDELQ $delin, "\n";
}

close (INDELQ);



exit;



