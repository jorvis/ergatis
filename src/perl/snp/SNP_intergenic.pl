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
my ($infile) = "GAK_contig.dir/gba_coords_com_name.doc";
my $motiffile;
my $pos;
my $help;

my $tag;

my $HELPTEXT = "Please read the README file located in this folder";

Getopt::Long::config("no_ignore_case");

my ($result) = GetOptions ("h|help"    => \$help,
			   "P=s"       => \$motiffile,
			   "T=s"       => \$tag);
			   

if ($result == 0) {
    die ("Command line parsing failed\n")
    }

if ($help) {
    print $HELPTEXT, "\n";
}

if (! defined $infile) {
    print STDERR "Cannot open $infile\n";
    exit(0);
}

if (! defined $motiffile) {
    print STDERR "Cannot open $pos\n";
    exit(0);
}

#### READ THE FILE INTO A ARRAY

my @infile = get_file_data($infile);
my @posarray = get_file_data($motiffile);


my $front = $tag . ".front";
my $ingene = $tag . ".ingene";
my $reject = $tag . ".reject";

#### CREATE HASHES OF THE ORF::END5::END3::COM_NAME FILE
#### ORF = KEYS 
#### VALUES = END5 or END3 or COM_NAME

######### OPEN FILE FOR PRINTOUT ##########

unless (open (FRONT, ">$front")){
    print STDERR "Cannot open output $front";
    exit(0);
}

print FRONT "LIST OF SNP LOCATED IN INTERGENIC REGION\n\n";

unless (open (IN, ">$ingene")) {
    print STDERR "Cannot open output $ingene";
    exit(0);
};


unless (open (REJ, ">$reject")) {
    print STDERR "Cannot open output $reject";
    exit(0);
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

    my @pos = split(" ", $posline);

    $pos = $pos[2];
    my $SNP_ID;
    ($SNP_ID = $pos[0]) =~ s/^>//;

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
	   
	    my $dist1 = $end5 - $end3m1;
	    my $dist2 = $end3p1 - $end3;
	    my $dist = $end5 - $pos;
	    print FRONT "SNP : $pos[0]  $pos[2] $pos[3] $pos[4] :\n";
	    print FRONT "This SNP is $dist bp in front of $orf\n";
	    print FRONT "|---$orfm1---> SNP $dist bp |---$orf---> <---$orfp1---|\n";
	    print FRONT "|---$orfm1---> $dist1 bp |---$orf---> $dist2 bp <---$orfp1---|\n";
	    print FRONT "$orfm1 ($lm1 bp): $ORF_COM{$orfm1}\n";
	    print FRONT "$orf ($lo bp) : $ORF_COM{$orf}\n";
	    print FRONT "$orfp1 ($lp1 bp): $ORF_COM{$orfp1}\n\n";
	     
	    last;

	      ### --> M -0-> -->

	  }elsif ($end5 < $end3 && $end5m1 < $end3m1 && $pos > $end3m1 && $end5p1 < $end3p1) {

	    print $pos, "--> S --> -->\n\n";
	   
	    my $dist1 = $end5 - $end3m1;
	    my $dist2 = $end5p1 - $end3;
	    my $dist = $end5 - $pos;
	    print FRONT "SNP : $pos[0]  $pos[2] $pos[3] $pos[4] :\n";
	    print FRONT "This SNP is $dist bp in front of $orf\n";
	    print FRONT "|---$orfm1---> SNP $dist bp |---$orf---> |---$orfp1--->\n";
	    print FRONT "|---$orfm1--->  $dist1 bp |---$orf---> $dist2 bp |---$orfp1--->\n";
	    print FRONT "$orfm1 ($lm1 bp): $ORF_COM{$orfm1}\n";
	    print FRONT "$orf ($lo bp) : $ORF_COM{$orf}\n";
	    print FRONT "$orfp1 ($lp1 bp): $ORF_COM{$orfp1}\n\n";

	    last;

	      ### <-- -M-> -0-> -->

	}elsif ($end5 < $end3 && $end5m1 < $end3m1 && $pos < $end3m1 && $pos > $end5m1 && $end5m2 > $end3m2 && $end5p1 < $end3p1) {

	    print $pos, "<--  -S-> --> -->\n\n";

	    my $dist1 = $end5 - $end3m1;
	    my $dist2 = $end5m1 - $end5m2;
	    my $dist3 = $end5p1 - $end3;
	    my $dist = $end5 - $pos;
	    my $dist4 = $end3m1 -$pos;
	    my $dist5 = ($pos - $len) - $end5m2;
	    print IN "SNP : $pos[0]  $pos[2] $pos[3] $pos[4] :\n";
	    print IN "This SNP is $dist bp in front of $orf\n";
	    print IN "but $dist4 bp from END3 within $orfm1\n";
	    print IN "and  $dist5 bp from END5 of $orfm2}\n";
	    print IN "<---$orfm2---| $dist5 bp |---$orfm1-SNP---> $dist bp |---$orf---> |---$orfp1--->\n";
	    print IN "<---$orfm2---| $dist2 bp |---$orfm1---> $dist1 bp |---$orf---> $dist3 bp|---$orfp1--->\n";
	    print IN "$orfm2 ($lm2 bp): $ORF_COM{$orfm2}\n";
	    print IN "$orfm1 ($lm1 bp): $ORF_COM{$orfm1}\n";
	    print IN "$orf ($lo bp) : $ORF_COM{$orf}\n";
	    print IN "$orfp1 ($lp1 bp): $ORF_COM{$orfp1}\n\n";
	   
	    last;

	       ### <-- -M-> -0-> <--
	    
	}elsif ($end5 < $end3 && $end5m1 < $end3m1 && $pos < $end3m1 && $pos > $end5m1 && $end5m2 > $end3m2 && $end5p1 > $end3p1) {

	    print $pos, "<--  -S-> --> <--\n\n";

	    my $dist1 = $end5 - $end3m1;
	    my $dist2 = $end5m1 - $end5m2;
	    my $dist3 = $end3p1 - $end3;
	    my $dist = $end5 - $pos;
	    my $dist4 = $end3m1 -$pos;
	    my $dist5 = ($pos - $len) - $end5m2;
	    print IN "SNP : $pos[0]  $pos[2] $pos[3] $pos[4] :\n";
	    print IN "This SNP is $dist bp in front of $orf\n";
	    print IN "but $dist4 bp from END3 within $orfm1\n";
	    print IN "and  $dist5 bp from END5 of $orfm2}\n";
	    print IN "<---$orfm2---| $dist5 bp |---$orfm1-SNP---> $dist bp |---$orf---> <---$orfp1---|\n";
	    print IN "<---$orfm2---| $dist2 bp |---$orfm1---> $dist1 bp |---$orf---> $dist3 bp <---$orfp1---|\n";
	    print IN "$orfm2 ($lm2 bp): $ORF_COM{$orfm2}\n";
	    print IN "$orfm1 ($lm1 bp): $ORF_COM{$orfm1}\n";
	    print IN "$orf ($lo bp) : $ORF_COM{$orf}\n";
	    print IN "$orfp1 ($lp1 bp): $ORF_COM{$orfp1}\n\n";
	   
	    last;

	       ### --> -M-> -0-> <--

	}elsif ($end5 < $end3 && $end5m1 < $end3m1 && $pos < $end3m1 && $pos > $end5m1 && $end5m2 < $end3m2 && $end5p1 > $end3p1) {

	    print $pos, "-->  -M-> --> <--\n\n";

	    my $dist1 = $end5 - $end3m1;
	    my $dist2 = $end5m1 - $end3m2;
	    my $dist3 = $end3p1 - $end3;
	    my $dist = $end5 - $pos;
	    my $dist4 = $end3m1 -$pos;
	    my $dist5 = $pos - $end5m2;
	    print IN "SNP : $pos[0]  $pos[2] $pos[3] $pos[4] :\n";
	    print IN "This SNP is $dist bp in front of $orf\n";
	    print IN "but $dist4 bp from END3 within $orfm1\n";
	    print IN "|---$orfm2---> |---$orfm1-SNP---> $dist bp |---$orf---> <---$orfp1---|\n";
	    print IN "|---$orfm2---> $dist2 bp |---$orfm1---> $dist1 bp |---$orf---> $dist3 bp <---$orfp1---|\n";
	    print IN "$orfm2 ($lm2 bp): $ORF_COM{$orfm2}\n";
	    print IN "$orfm1 ($lm1 bp): $ORF_COM{$orfm1}\n";
	    print IN "$orf ($lo bp) : $ORF_COM{$orf}\n";
	    print IN "$orfp1 ($lp1 bp): $ORF_COM{$orfp1}\n\n";
	    
	    last;

	       ### --> -M-> -0-> -->

	}elsif ($end5 < $end3 && $end5m1 < $end3m1 && $pos < $end3m1 && $pos > $end5m1 && $end5m2 < $end3m2 && $end5p1 < $end3p1) {

	    print $pos, "-->  -M-> --> -->\n\n";

	    my $dist1 = $end5 - $end3m1;
	    my $dist2 = $end5m1 - $end3m2;
	    my $dist3 = $end5p1 - $end3;
	    my $dist = $end5 - $pos;
	    my $dist4 = $end3m1 -$pos;
	    print IN "MOTIF : $pos[1] at position $pos[0]:\n";
	    print IN "This motif is $dist bp in front of $orf\n";
	    print IN "but $dist4 bp from END3 within $orfm1\n";
	    print IN "|---$orfm2---> |---$orfm1-MOTIF---> $dist bp |---$orf---> |---$orfp1--->\n";
	    print IN "|---$orfm2---> $dist2 bp |---$orfm1---> $dist1 bp |---$orf---> $dist3 bp |---$orfp1--->\n";
	    print IN "$orfm2 ($lm2 bp): $ORF_COM{$orfm2}\n";
	    print IN "$orfm1 ($lm1 bp): $ORF_COM{$orfm1}\n";
	    print IN "$orf ($lo bp) : $ORF_COM{$orf}\n";
	    print IN "$orfp1 ($lp1 bp): $ORF_COM{$orfp1}\n\n";
	    
	    last;


	    #### GENE -1 IN REVERSE ORIENTATION


	       ### <-- <-- M -0-> -->

	}elsif ($end5 < $end3 && $end5m1 > $end3m1 && $pos > $end5m1 && $end5m2 > $end3m2 && $end5p1 < $end3p1) {

	     print $pos, "<--  <-- M  --> -->\n\n";

	     my $dist1 = $end3m1 -$end5m2;
	     my $dist2 = $end5 - $end5m1;
	     my $dist3 = $end5p1 - $end3;
	     my $dist = $end5 - $pos;
	     my $dist4 = ($pos - $len) - $end5m1;
	     print FRONT "SNP : $pos[0]  $pos[2] $pos[3] $pos[4] :\n";
	     print FRONT "This SNP is $dist bp in front of $orf\n";
	     print FRONT "and $dist4 from END5 of $orfm1\n";
	     print FRONT "<---$orfm2---|<---$orfm1---| $dist4 bp SNP $dist bp |---$orf---> |---$orfp1--->\n";
	     print FRONT "<---$orfm2---| $dist1 bp <---$orfm1---| $dist2 bp |---$orf---> $dist3 bp |---$orfp1--->\n";
	     print FRONT "$orfm2 ($lm2 bp): $ORF_COM{$orfm2}\n";
	     print FRONT "$orfm1 ($lm1 bp): $ORF_COM{$orfm1}\n";
	     print FRONT "$orf ($lo bp) : $ORF_COM{$orf}\n";
	     print FRONT "$orfp1 ($lp1 bp): $ORF_COM{$orfp1}\n\n";

	     last;

	        ### <-- <-- M --> <--

	 }elsif ($end5 < $end3 && $end5m1 > $end3m1 && $pos > $end5m1 && $end5m2 > $end3m2 && $end5p1 > $end3p1) {

	     print $pos, "<--  <-- M  --> <--\n\n";

	     my $dist1 = $end3m1 -$end5m2;
	     my $dist2 = $end5 - $end5m1;
	     my $dist3 = $end3p1 - $end3;
	     my $dist = $end5 - $pos;
	     my $dist4 = ($pos - $len) - $end5m1;
	     print FRONT "SNP : $pos[0]  $pos[2] $pos[3] $pos[4] :\n";
	     print FRONT "This SNP is $dist bp in front of $orf\n";
	     print FRONT "and $dist4 from END5 of $orfm1\n";
	     print FRONT "<---$orfm2---|<---$orfm1---| $dist4 bp SNP $dist bp |---$orf---> <---$orfp1---|\n";
	     print FRONT "<---$orfm2---| $dist1 bp <---$orfm1---| $dist2 bp |---$orf---> $dist3 bp <---$orfp1---|\n";
	     print FRONT "$orfm2 ($lm2 bp): $ORF_COM{$orfm2}\n";
	     print FRONT "$orfm1 ($lm1 bp): $ORF_COM{$orfm1}\n";
	     print FRONT "$orf ($lo bp) : $ORF_COM{$orf}\n";
	     print FRONT "$orfp1 ($lp1 bp): $ORF_COM{$orfp1}\n\n";

	     last;

	        ### --> <-- M --> <--

	}elsif ($end5 < $end3 && $end5m1 > $end3m1 && $pos > $end5m1 && $end5m2 < $end3m2 && $end5p1 > $end3p1) {

	     print $pos, "-->  <-- M  --> <--\n\n";

	     my $dist1 = $end3m1 -$end3m2;
	     my $dist2 = $end5 - $end5m1;
	     my $dist3 = $end3p1 - $end3;
	     my $dist = $end5 - $pos;
	     my $dist4 = ($pos - $len) - $end5m1;
	     print FRONT "SNP : $pos[0]  $pos[2] $pos[3] $pos[4] :\n";
	     print FRONT "This SNP is $dist bp in front of $orf\n";
	     print FRONT "and $dist4 from END5 of $orfm1\n";
	     print FRONT "|---$orfm2---><---$orfm1---| $dist4 bp SNP $dist bp |---$orf---> <---$orfp1---|\n";
	     print FRONT "|---$orfm2---> $dist1 bp <---$orfm1---| $dist2 bp |---$orf---> $dist3 bp <---$orfp1---|\n";
	     print FRONT "$orfm2 ($lm2 bp): $ORF_COM{$orfm2}\n";
	     print FRONT "$orfm1 ($lm1 bp): $ORF_COM{$orfm1}\n";
	     print FRONT "$orf ($lo bp) : $ORF_COM{$orf}\n";
	     print FRONT "$orfp1 ($lp1 bp): $ORF_COM{$orfp1}\n\n";

	     last;

	        ### --> <-- M --> -->

	 }elsif ($end5 < $end3 && $end5m1 > $end3m1 && $pos > $end5m1 && $end5m2 < $end3m2 && $end5p1 < $end3p1) {

	     print $pos, "-->  <-- M  --> -->\n\n";

	     my $dist1 = $end3m1 -$end3m2;
	     my $dist2 = $end5 - $end5m1;
	     my $dist3 = $end5p1 - $end3;
	     my $dist = $end5 - $pos;
	     my $dist4 = ($pos - $len) - $end5m1;
	     print FRONT "SNP : $pos[0]  $pos[2] $pos[3] $pos[4] :\n";
	     print FRONT "This SNP is $dist bp in front of $orf\n";
	     print FRONT "and $dist4 from END5 of $orfm1\n";
	     print FRONT "|---$orfm2---><---$orfm1---| $dist4 bp SNP $dist bp |---$orf---> |---$orfp1--->\n";
	     print FRONT "|---$orfm2---> $dist1 bp <---$orfm1---| $dist2 bp |---$orf---> $dist3 bp |---$orfp1--->\n";
	     print FRONT "$orfm2 ($lm2 bp): $ORF_COM{$orfm2}\n";
	     print FRONT "$orfm1 ($lm1 bp): $ORF_COM{$orfm1}\n";
	     print FRONT "$orf ($lo bp) : $ORF_COM{$orf}\n";
	     print FRONT "$orfp1 ($lp1 bp): $ORF_COM{$orfp1}\n\n";

	     last;  


	### B.  END5 > END3 GENE 0 IN REVERSE ORIENTATION

	        ### --> <-M- <-

	}elsif ($end5 > $end3 && $end5m1 < $end3m1 && $pos > $end3 && $end5p1 > $end3p1) {
	
	    print $pos, "--> <-M- <--\n\n";

	    
	    my $dist1 = $end3 - $end3m1;
	    my $dist2 = $end3p1 - $end5;
	    print REJ "MOTIF : $pos[1] at position $pos[0]:\n";
	    print REJ "This motif is in $orf\n";
	    print REJ "|---$orfm1---> <---$orf-MOTIF---| <---$orfp1---|\n";
	    print REJ "|---$orfm1---> $dist1 bp <---$orf-MOTIF---| $dist2 bp <---$orfp1---|\n";
	    print REJ "$orfm1 ($lm1 bp): $ORF_COM{$orfm1}\n";
	    print REJ "$orf ($lo bp) : $ORF_COM{$orf}\n";
	    print REJ "$orfp1 ($lp1 bp): $ORF_COM{$orfp1}\n\n";

	   last;
 
	       ### --> <-M- -->

	}elsif ($end5 > $end3 && $end5m1 < $end3m1 && $pos > $end3 && $end5p1 < $end3p1) {
	
	    print $pos, "--> <-M- -->\n\n";

	    my $dist3 = $end5p1 - $pos;;
	    my $dist1 = $end3 - $end3m1;
	    my $dist2 = $end5p1 - $end5;
	    print IN "MOTIF : $pos[1] at position $pos[0]:\n";
	    print IN "This motif is in $orf\n";
	    print IN "|---$orfm1--->  <---$orf-MOTIF---| $dist3 bp |---$orfp1--->\n";
	    print IN "|---$orfm1---> $dist1 bp <---$orf-MOTIF---| $dist2 bp |---$orfp1--->\n";
	    print IN "$orfm1 ($lm1 bp): $ORF_COM{$orfm1}\n";
	    print IN "$orf ($lo bp) : $ORF_COM{$orf}\n";
	    print IN "$orfp1 ($lp1 bp): $ORF_COM{$orfp1}\n\n";

	    last;	    

	       ### --> M <-0-

	}elsif ($end5 > $end3 && $end5m1 < $end3m1 && $pos < $end3 && $pos > $end3m1) {
	
	    print $pos, "--> M <--\n\n";

	    my $dist = $end3 - $pos;
	    my $dist2 = ($pos -$len) - $end3m1;
	    my $dist3 = $end3 - $end3m1;
	    print REJ "MOTIF : $pos[1] at position $pos[0]:\n";
	    print REJ "This motif is in between the END3 of $orfm1 and $orf\n";
	    print REJ "|---$orfm1---> $dist2 bp MOTIF $dist bp <---$orf---|\n";
	    print REJ "|---$orfm1---> $dist3 bp <---$orf---|\n";
	    print REJ "$orfm1 ($lm1 bp): $ORF_COM{$orfm1}\n";
	    print REJ "$orf ($lo bp) : $ORF_COM{$orf}\n\n";

	    last;    

	       ### <-- -M-> <-0-

	}elsif ($end5 > $end3 && $end5m1 < $end3m1 && $pos < $end3m1 && $end5m2 > $end3m2) {

	     print $pos, "<-- -M-> <--\n\n";

	     my $dist = ($pos - $len) - $end5m2;
	     my $dist1 = $end3 - $end3m1;
	     my $dist2 = $end5m1 - $end5m2; 
	     print IN "MOTIF : $pos[1] at position $pos[0]:\n";
	     print IN "This motif is in $orfm1\n";
	     print IN "and $dist bp from $orfm2\n";
	     print IN "<---$orfm2---| $dist bp |---$orfm1-MOTIF---> <---$orf---|\n";
	     print IN "<---$orfm2---| $dist2 bp |---$orfm1-MOTIF---> $dist1 bp <---$orf---|\n";
	     print IN "$orfm2 ($lm2 bp): $ORF_COM{$orfm2}\n";
	     print IN "$orfm1 ($lm1 bp): $ORF_COM{$orfm1}\n";
	     print IN "$orf ($lo bp) : $ORF_COM{$orf}\n\n";

	     last; 

	        ### --> -M-> <-0-

	 }elsif ($end5 > $end3 && $end5m1 < $end3m1 && $pos < $end3m1 && $end5m2 < $end3m2) {

	     print $pos, "--> -M-> <--\n\n";

	     my $dist1 = $end5m1 - $end3m2;
	     my $dist2 = $end3 - $end3m1; 
	     print REJ "MOTIF : $pos[1] at position $pos[0]:\n";
	     print REJ "This motif is in $orfm1\n";
	     print REJ "|---$orfm2---> |---$orfm1-MOTIF---> <---$orf---|\n";
	     print REJ "|---$orfm2---> $dist1 bp |---$orfm1-MOTIF---> $dist2 bp <---$orf---|\n";
	     print REJ "$orfm2 ($lm2 bp): $ORF_COM{$orfm2}\n";
	     print REJ "$orfm1 ($lm1 bp): $ORF_COM{$orfm1}\n";
	     print REJ "$orf ($lo bp) : $ORF_COM{$orf}\n\n";


	     last; 

	        ### <-- <-M0- <--

	 }elsif ($end5 > $end3 && $end5m1 > $end3m1 && $end5p1 > $end3p1 && $pos > $end3) {

	    print $pos, "<-- <-M- <--\n\n";

	    my $dist = ($pos - $len) - $end5m1;
	    my $dist1 = $end3 - $end5m1;
	    my $dist2 = $end3p1 - $end5;
	    print IN "MOTIF : $pos[1] at position $pos[0]:\n";
	    print IN "This motif is in $orf\n";
	    print IN "<---$orfm1---| $dist bp <---$orf-MOTIF---| <---$orfp1---|\n";
	    print IN "<---$orfm1---| $dist1 bp <---$orf-MOTIF---| $dist2 bp <---$orfp1---|\n";
	    print IN "$orfm1 ($lm1 bp): $ORF_COM{$orfm1}\n";
	    print IN "$orf ($lo bp) : $ORF_COM{$orf}\n";
	    print IN "$orfp1 ($lp1 bp): $ORF_COM{$orfp1}\n\n";

	    last;

	       ### <-- <-M0- -->

	}elsif ($end5 > $end3 && $end5m1 > $end3m1 && $end5p1 < $end3p1 && $pos > $end3) {

	    print $pos, "<-- <-M- -->\n\n";

	    my $dist = $end5p1 - $pos;
	    my $dist3 = ($pos - $len) - $end5m1;
	    my $dist1 = $end3 - $end5m1;
	    my $dist2 = $end5p1 - $end5;
	    print IN "MOTIF : $pos[1] at position $pos[0]:\n";
	    print IN "This motif is in $orf\n";
	    print IN "<---$orfm1---| $dist3 bp <---$orf-MOTIF---| $dist bp |---$orfp1--->\n";
	    print IN "<---$orfm1---| $dist1 bp <---$orf-MOTIF---| $dist2 bp |---$orfp1--->\n";
	    print IN "$orfm1 ($lm1 bp): $ORF_COM{$orfm1}\n";
	    print IN "$orf ($lo bp) : $ORF_COM{$orf}\n";
	    print IN "$orfp1 ($lp1 bp): $ORF_COM{$orfp1}\n\n";


	    last;

	       ### <-- <-- M <-0-

	}elsif ($end5 > $end3 && $end5m1 > $end3m1 && $pos < $end3 && $end5m2 > $end3m2) {

	    print $pos, "<-- <-- M <--\n\n";

	    my $dist = ($pos - $len) - $end5m1;
	    my $dist1 = $end3 - $end5m1;
	    my $dist2 = $end3m1 - $end5m2;
	    print FRONT "SNP : $pos[0]  $pos[2] $pos[3] $pos[4] :\n";
	    print FRONT "This SNP is $dist bp from END5 of $orfm1\n";
	    print FRONT "<---$orfm2---| <---$orfm1--| $dist bp SNP <---$orf---|\n";
	    print FRONT "<---$orfm2---| $dist 2 bp <---$orfm1--| $dist1 bp <---$orf---|\n";
	    print FRONT "$orfm2 ($lm2 bp): $ORF_COM{$orfm2}\n";
	    print FRONT "$orfm1 ($lm1 bp): $ORF_COM{$orfm1}\n";
	    print FRONT "$orf ($lo bp) : $ORF_COM{$orf}\n\n";
      
	    last;

	       ### --> <-- M <-0-

	}elsif ($end5 > $end3 && $end5m1 > $end3m1 && $pos < $end3 && $end5m2 < $end3m2) {

	    print $pos, "--> <-- M <--\n\n";

	    my $dist = ($pos - $len) - $end5m1;
	    my $dist1 = $end3 - $end5m1;
	    my $dist2 = $end3m1 - $end3m2;
	    print FRONT "SNP : $pos[0]  $pos[2] $pos[3] $pos[4] :\n";
	    print FRONT "This SNP is $dist bp from END5 of $orfm1\n";
	    print FRONT "|---$orfm2---> <---$orfm1--| $dist bp SNP <---$orf---|\n";
	    print FRONT "|---$orfm2---> $dist 2 bp <---$orfm1--| $dist1 bp <---$orf---|\n";
	    print FRONT "$orfm2 ($lm2 bp): $ORF_COM{$orfm2}\n";
	    print FRONT "$orfm1 ($lm1 bp): $ORF_COM{$orfm1}\n";
	    print FRONT "$orf ($lo bp) : $ORF_COM{$orf}\n\n";
      
	    last;


	}
    }
	    
$i++

}  ### END OF FOR EACH $e5 < END5

}  ### END OF FOREACH GOING THROUGH LIST OF MOTIFS POSITION


close (FRONT);
close (IN);
close (REJ);


system("mv $tag.front  SNP_intergenic.txt");

#system("rm $tag.front");
#system("rm $tag.ingene");
#system("rm $tag.reject");


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







	





















