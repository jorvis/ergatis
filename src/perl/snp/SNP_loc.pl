#! /usr/local/bin/perl -w


###############################################################
### MOTIF_POSITION_FINDER FOR INTERGENIC SNP
### GIVEN A POSITION FIND THE GENE IT IS LOCATED IN FRONT OF
################################################################

use strict;
use Getopt::Long;


my (%end5ORF);
my (%end3ORF);
my (%ORF_COM);
my ($infile); ## coords file
my $motiffile; ## list of position (output of affy2.pl)
my $pos;
my $newline;
my @SNP_to_gene;
my $refpos;

my $help;

my $tag;

my $HELPTEXT = " 
USAGE: SNP_location.pl -F SNP_list -C coords_com_name file 
 
SNP_list: List of SNP (output from affy2.pl)

coords_com_name: Coordinate file with com_name separated by ::

-h Help file 

For for Info please E-mail jravel\@tigr.org

";

Getopt::Long::config("no_ignore_case");

my ($result) = GetOptions ("h|help"    => \$help,
			   "F=s"       => \$motiffile,
			   "C=s"       => \$infile,
			   "P=s"       => \$refpos,
			   "O=s"       => \$tag);
			   


if ($result == 0) {
    die ("Command line parsing failed\n")
    }

if ($help) {
    print $HELPTEXT, "\n";
    exit(0);
}

if (! defined $infile) {
    print STDERR "Cannot open $infile\n";
    exit(0);
}


if (! defined $motiffile) {
    print STDERR "Cannot open $motiffile\n";
    exit(0);
}


$infile = "/home/jravel/src/SNP/" . $infile;



#### READ THE FILE INTO A ARRAY

my @infile = get_file_data($infile);
my @posarray = get_file_data($motiffile);


my $outputfile = $tag . "_location.txt";


#### CREATE HASHES OF THE ORF::END5::END3::COM_NAME FILE
#### ORF = KEYS 
#### VALUES = END5 or END3 or COM_NAME

if ($tag eq "SNP") {
my $header = "SNP_id\tREF LENGTH\tREF_POS\tSNP\tQUE_POS\tORF_ID\tCOM_NAME";
push(@SNP_to_gene, $header);
}elsif ($tag eq "INDEL") {

my $header = "INDEL_id\tREF POS\tINDEL_LEN\tQUE_POS\tQUE_ID\tORF_ID\tCOM_NAME";
push(@SNP_to_gene, $header);
}


foreach my $line (@infile) { ## READ THE GENE COORDS FILE 

    chomp $line;

    ## SPLIT EACH FIELD

    my @field = split("::", $line);

    ## IF NO COM_NAME, COM_NAME = BLANK

    if (! defined $field[3]) {
	$field[3] = " ";
    }

    my $baorf = $field[0];
    my $end5 = $field[1];
    my $end3 = $field[2];
    (my $com_name = $field[3]) =~ s/,/ /g;
    

    ### CREATE THE HASHES
    
    $ORF_COM{$baorf} = $com_name;  ### KEYS = ORF_NAME VALUES = COM_NAME
    $end5ORF{$baorf} = $end5;   ### KEYS = ORF_NAME VALUES = END5
    $end3ORF{$baorf} = $end3;  ### KEYS = ORF_NAME VALUES = END3


} #end of foreach statement



## SORT THE HASH NUMERICALLY ON THE END5

my @end5 = sort { $a <=> $b} values %end5ORF;


## REVERSE THE HASH (KEYS = END5 VALUES=ORFs)

my %R_end5ORF = reverse %end5ORF;







### GO THROUGH THE POSITION FILE - TAKE EACH LINE (SPLIT) AND USE POSITION = $pos 
### INTO THE NEXT FOREACH STATEMENT (GO THROUGH THE HASH OF ORF-ENDs)

foreach my $posline (@posarray) { ### LIST OF SNP POSITION ON GBA

    chomp $posline;

    my @pos = split(' ', $posline);

    $pos = $pos[$refpos];

    my  $len = 1; ### FOR SNP THIS IS ONE


my $i = 0;

foreach my $e5 (@end5){


### THE MAJOR TEST IS FOR THE POSITION OF THE SNP TO BE UPSTREAM FROM AN END5 
### THERE ARE TWO POSSIBILITY THE GENE (WHICH THE END5 BELONGS TO) IS IN FORWARD ORIENTATION OR REVERSE ORIENTATION
### IF IN REVERSE - IS THE MOTIF IN THE MIDDLE OF THE GENE - BUT IS THE MOTIF IN FRONT OF THE NEXT GENE (FORWARD/REVERSE)
### IF IN FORWARD - IS THE MOTIF IN THE PREVIOUS GENE (FORWARD OR REVERSE)
### THIS SCRIPT TAKES ACCOUNT OF ALL THE POSSIBILITIES - AND CALCULATES THE DISTANCE FROM THE MOTIF TO END5s AND
### GIVES THE NAME OF THE GENE ITS IN FRONT OF

 
    if ($pos < $e5) {  ## IF POS < END5  



my $end5 = $e5;
my $end3 = $end3ORF{$R_end5ORF{$e5}};
my $end5p1 = $end5[$i+1];
my $end3p1 = $end3ORF{$R_end5ORF{$end5[$i+1]}};
my $end5m1 = $end5[$i-1];
my $end3m1 = $end3ORF{$R_end5ORF{$end5[$i-1]}};
my $end5m2 = $end5[$i-2];
my $end3m2 = $end3ORF{$R_end5ORF{$end5[$i-2]}};
my $orf = $R_end5ORF{$e5};
my $orfp1 = $R_end5ORF{$end5[$i+1]};
my $orfm1 = $R_end5ORF{$end5[$i-1]};
my $orfm2 = $R_end5ORF{$end5[$i-2]};
my $lm2 = abs($end5m2-$end3m2);
my $lm1 = abs($end5m1-$end3m1);
my $lo = abs($end5-$end3);
my $lp1 = abs($end5p1-$end3p1);
my $com_m2 =  $ORF_COM{$orfm2};
my $com_m1 =  $ORF_COM{$orfm1};
my $com =  $ORF_COM{$orf};
my $com_p1 = $ORF_COM{$orfp1};
	
#### A. END5 < END3  GENE 0 IN FORWARD ORIENTATION  -->


	### A.1  END5m1 < END3m1 GENE -1 IN FORWARD ORIENTATION --> -O-> 

	     ### --> M -0-> <--

	if ($end5 < $end3 && $end5m1 < $end3m1 && $pos > $end3m1 && $end5p1 > $end3p1) {

	    print $pos, "--> S --> <--\n\n";
	    
	    $newline = "$pos[0]\t$pos[1]\t$pos[2]\t$pos[3]\t$pos[4]\tIntergenic";
	    push(@SNP_to_gene, $newline);
	     
	    last;

	      ### --> M -0-> -->

	  }elsif ($end5 < $end3 && $end5m1 < $end3m1 && $pos > $end3m1 && $end5p1 < $end3p1) {

	    print $pos, "--> S --> -->\n\n";

	    $newline = "$pos[0]\t$pos[1]\t$pos[2]\t$pos[3]\t$pos[4]\tIntergenic";
	    push(@SNP_to_gene, $newline);

	    last;

	      ### <-- -M-> -0-> -->

	}elsif ($end5 < $end3 && $end5m1 < $end3m1 && $pos < $end3m1 && $pos > $end5m1 && $end5m2 > $end3m2 && $end5p1 < $end3p1) {

	    print $pos, "<--  -S-> --> -->\n\n";


	    $newline = "$pos[0]\t$pos[1]\t$pos[2]\t$pos[3]\t$pos[4]\t" . $orfm1. "\t".  $ORF_COM{$orfm1};
	    push(@SNP_to_gene, $newline);
	   
	    last;

	       ### <-- -M-> -0-> <--
	    
	}elsif ($end5 < $end3 && $end5m1 < $end3m1 && $pos < $end3m1 && $pos > $end5m1 && $end5m2 > $end3m2 && $end5p1 > $end3p1) {

	    print $pos, "<--  -S-> --> <--\n\n";

	   
	    $newline = "$pos[0]\t$pos[1]\t$pos[2]\t$pos[3]\t$pos[4]\t" . $orfm1. "\t".  $ORF_COM{$orfm1};
	    push(@SNP_to_gene, $newline);

	    last;

	       ### --> -M-> -0-> <--

	}elsif ($end5 < $end3 && $end5m1 < $end3m1 && $pos < $end3m1 && $pos > $end5m1 && $end5m2 < $end3m2 && $end5p1 > $end3p1) {

	    print $pos, "-->  -M-> --> <--\n\n";

	    
	    $newline = "$pos[0]\t$pos[1]\t$pos[2]\t$pos[3]\t$pos[4]\t" . $orfm1. "\t".  $ORF_COM{$orfm1};
	    push(@SNP_to_gene, $newline);


	    last;

	       ### --> -M-> -0-> -->

	}elsif ($end5 < $end3 && $end5m1 < $end3m1 && $pos < $end3m1 && $pos > $end5m1 && $end5m2 < $end3m2 && $end5p1 < $end3p1) {

	    print $pos, "-->  -M-> --> -->\n\n";

	    
	    $newline = "$pos[0]\t$pos[1]\t$pos[2]\t$pos[3]\t$pos[4]\t" . $orfm1. "\t".  $ORF_COM{$orfm1};
	    push(@SNP_to_gene, $newline);

	    last;


	    #### GENE -1 IN REVERSE ORIENTATION


	       ### <-- <-- M -0-> -->

	}elsif ($end5 < $end3 && $end5m1 > $end3m1 && $pos > $end5m1 && $end5m2 > $end3m2 && $end5p1 < $end3p1) {

	     print $pos, "<--  <-- M  --> -->\n\n";


	     $newline = "$pos[0]\t$pos[1]\t$pos[2]\t$pos[3]\t$pos[4]\tIntergenic";
	     push(@SNP_to_gene, $newline);

	     last;

	        ### <-- <-- M --> <--

	 }elsif ($end5 < $end3 && $end5m1 > $end3m1 && $pos > $end5m1 && $end5m2 > $end3m2 && $end5p1 > $end3p1) {

	     print $pos, "<--  <-- M  --> <--\n\n";


	     $newline = "$pos[0]\t$pos[1]\t$pos[2]\t$pos[3]\t$pos[4]\tIntergenic";
	     push(@SNP_to_gene, $newline);

	     last;

	        ### --> <-- M --> <--

	}elsif ($end5 < $end3 && $end5m1 > $end3m1 && $pos > $end5m1 && $end5m2 < $end3m2 && $end5p1 > $end3p1) {

	     print $pos, "-->  <-- M  --> <--\n\n";


	     $newline = "$pos[0]\t$pos[1]\t$pos[2]\t$pos[3]\t$pos[4]\tIntergenic";
	     push(@SNP_to_gene, $newline);

	     last;

	        ### --> <-- M --> -->

	 }elsif ($end5 < $end3 && $end5m1 > $end3m1 && $pos > $end5m1 && $end5m2 < $end3m2 && $end5p1 < $end3p1) {

	     print $pos, "-->  <-- M  --> -->\n\n";


	     $newline = "$pos[0]\t$pos[1]\t$pos[2]\t$pos[3]\t$pos[4]\tIntergenic";
	     push(@SNP_to_gene, $newline);

	     last;  


	### B.  END5 > END3 GENE 0 IN REVERSE ORIENTATION

	        ### --> <-M- <-

	}elsif ($end5 > $end3 && $end5m1 < $end3m1 && $pos > $end3 && $end5p1 > $end3p1) {
	
	    print $pos, "--> <-M- <--\n\n";

	    $newline = "$pos[0]\t$pos[1]\t$pos[2]\t$pos[3]\t$pos[4]\t" . $orf. "\t".  $ORF_COM{$orf};
	    push(@SNP_to_gene, $newline);

	    last;
 
	       ### --> <-M- -->

	}elsif ($end5 > $end3 && $end5m1 < $end3m1 && $pos > $end3 && $end5p1 < $end3p1) {
	
	    print $pos, "--> <-M- -->\n\n";


	    $newline = "$pos[0]\t$pos[1]\t$pos[2]\t$pos[3]\t$pos[4]\t" . $orf. "\t".  $ORF_COM{$orf};
	    push(@SNP_to_gene, $newline);

	    last;	    

	       ### --> M <-0-

	}elsif ($end5 > $end3 && $end5m1 < $end3m1 && $pos < $end3 && $pos > $end3m1) {
	
	    print $pos, "--> M <--\n\n";


	    $newline = "$pos[0]\t$pos[1]\t$pos[2]\t$pos[3]\t$pos[4]\tIntergenic";
	    push(@SNP_to_gene, $newline);

	    last;    

	       ### <-- -M-> <-0-

	}elsif ($end5 > $end3 && $end5m1 < $end3m1 && $pos < $end3m1 && $end5m2 > $end3m2) {

	     print $pos, "<-- -M-> <--\n\n";

	     $newline = "$pos[0]\t$pos[1]\t$pos[2]\t$pos[3]\t$pos[4]\t" . $orfm1. "\t".  $ORF_COM{$orfm1};
	     push(@SNP_to_gene, $newline);

	     last; 

	        ### --> -M-> <-0-

	 }elsif ($end5 > $end3 && $end5m1 < $end3m1 && $pos < $end3m1 && $end5m2 < $end3m2) {

	     print $pos, "--> -M-> <--\n\n";


	     $newline = "$pos[0]\t$pos[1]\t$pos[2]\t$pos[3]\t$pos[4]\t" . $orfm1. "\t".  $ORF_COM{$orfm1};
	     push(@SNP_to_gene, $newline);


	     last; 

	        ### <-- <-M0- <--

	 }elsif ($end5 > $end3 && $end5m1 > $end3m1 && $end5p1 > $end3p1 && $pos > $end3) {

	    print $pos, "<-- <-M- <--\n\n";

	    $newline = "$pos[0]\t$pos[1]\t$pos[2]\t$pos[3]\t$pos[4]\t" . $orf. "\t".  $ORF_COM{$orf};
	    push(@SNP_to_gene, $newline);

	    last;

	       ### <-- <-M0- -->

	}elsif ($end5 > $end3 && $end5m1 > $end3m1 && $end5p1 < $end3p1 && $pos > $end3) {

	    print $pos, "<-- <-M- -->\n\n";

	    $newline = "$pos[0]\t$pos[1]\t$pos[2]\t$pos[3]\t$pos[4]\t" . $orf. "\t".  $ORF_COM{$orf};
	    push(@SNP_to_gene, $newline);


	    last;

	       ### <-- <-- M <-0-

	}elsif ($end5 > $end3 && $end5m1 > $end3m1 && $pos < $end3 && $end5m2 > $end3m2) {

	    print $pos, "<-- <-- M <--\n\n";

	    $newline = "$pos[0]\t$pos[1]\t$pos[2]\t$pos[3]\t$pos[4]\tIntergenic";
	    push(@SNP_to_gene, $newline);

	    last;

	       ### --> <-- M <-0-

	}elsif ($end5 > $end3 && $end5m1 > $end3m1 && $pos < $end3 && $end5m2 < $end3m2) {

	    print $pos, "--> <-- M <--\n\n";
      
	    $newline = "$pos[0]\t$pos[1]\t$pos[2]\t$pos[3]\t$pos[4]\tIntergenic";
	    push(@SNP_to_gene, $newline);


	    last;


	}
    }
	    
$i++

}  ### END OF FOR EACH $e5 < END5

}  ### END OF FOREACH GOING THROUGH LIST OF MOTIFS POSITION

;

unless (open(FILEOUT, ">$outputfile")) {
    print STDERR "Can't open file!!\n\n";
    exit;
}

foreach my $line2 (@SNP_to_gene) {
    print FILEOUT $line2, "\n";
}

close (FILEOUT);



exit;

################################################################
##
## SUBROUTINE  get_file_data
##
## A subroutine to get data from a file given its filename
################################################################


sub get_file_data {

    my($filename) = @_;

    use strict;
    use warnings;

    # Initialize variables
    my @filedata = (  );

    unless( open(GET_FILE_DATA, $filename) ) {
        print STDERR "Cannot open file \"$filename\"\n\n";
        exit;
    }

    @filedata = <GET_FILE_DATA>;

    close GET_FILE_DATA;

    return @filedata;
}







	





















