#! /usr/local/bin/perl -w

############################################################################
###
###                       index_get_qual.pl
###
###   This script retrieve information from a .tcov file using an index 
###   created with get_index_qual.pl
###
############################################################################


use strict;
use lib "/home/jravel/lib/";
use sub;
use Statistics::Descriptive;
use Getopt::Long;
$| = 1;

my $end;
my $pos;
my $help;
my $newcov;
my @quality;
my @qualityQ;
my $tcov2 = "/usr/local/projects/anthrax/jravel/GBA_cov.dir/ungapped.tcov";
my $tcov = "/usr/local/projects/anthrax/jravel/index.tcov";
my $list;
my @coverage;
my $fh;
my @list;
my $fh2;
my $record;
my $record2;
my $record3;
my $indlen;
my $end5R;
my $end3R;
my $end5Q;
my $end3Q;
my $length;
my @cov3;
my $tag;
my @coverage3;
my $lencov;
my $covref;
my $qualref;
my $lenque;
my $covque;
my $qualque;
my $refscore;
my $quescore;
my $insscore;
my $refcov;
my $quecov;
my $median;
my $medianQ;
my %locat;
my @location;
my @loca;
my $indel_id;

Getopt::Long::config("no_ignore_case");

my $result = GetOptions ("h|help"     => \$help,
			 "T=s"        => \$tag,
			 "F=s"        => \$list);

my $HELPTEXT = qq~
USAGE -F list file (INDEL_info.txt)

The list file is the output of the SNP pipeline and looks like this:
This is the output of indel_info.pl 


INDEL-REF-1:626874 1 bp ID94 (123 bp) INS on REF
R: 626854 - 626895 Q: 845 - 885
R: CTAGTAATGGAATTAATGCAACGTGGAACATGTAATGTAACA
Q: CTAGTAATGGAATTAATGCA-CGTGGAACATGTAATGTAACA
\/\/
INDEL-REF-2:748874 1 bp ID122 (123 bp) INS on REF
R: 748854 - 748895 Q: 970 - 1010
R: AGTATTATGTTAGATTTCGAATATAACAATATTTATTATTAA
Q: GTATTATGTTAGATTTTCGA-TATAACAATATTTATTATAAT
\/\/
...

INDEL-QUE-1:102159 1 bp ID121revcom (123 bp) INS on QUE
R: 102139 - 102179 Q: 1635 - 1676
R: TACAATCGATGAATTAAAAG-AACGTGGACTTTGGATTGCTG
Q: ACAGAAAGGAACTATCTTATCATTCACGAATTAATTTAATCC
\/\/
INDEL-QUE-2:137766 1 bp ID33revcom (123 bp) INS on QUE
R: 137746 - 137786 Q: 243 - 284
R: ACATGAAAAAGGATTAGAGC-AAGGGATATACATTTCTATAC
Q: ACATGAAAAAGGATTAGAGCCAAGGGATATACATTTCTATAC
\/\/



The offset (index.tcov file )is used in this script to retrieve information from the .tcov file.
Required is the file INDEL_location.txt (from the pipeline).
Perl function tell is used to grab the offset of each line. The seek function in index_get_qual.pl
is used to position the pointer to a particular line for retrieval in a large file.


HELP: -h (This text)
For questions and report bug please Email : jravel\@tigr.org


~;


if ($result == 0) {
    die ("Command line parsing failed\n")
    }

if ($help) {
    print $HELPTEXT, "\n";
    exit;
}


if (! defined $list) {
    print STDERR "You must specify a file name (-F)\n";
    exit(0);
}

if (! defined $tag) {
    print STDERR "You must specify a folder name (-T)\n";
    exit(0);
}


##### DO THE TRANSFORMATION FROM 6615 to 6609 ACCORDING TO:

############################################################################
# asmbl_id    asm_lend    asm_rend    sub_asmbl_id sub_asm_lend sub_asm_rend
# ----------- ----------- ----------- ------------ ------------ ------------
#       6615           1     4651062         6609       576236      5227297
#       6615     4651063     5227297         6609            1       576235
#
#       INCLUDING THESE LOCAL DIFFERENCES
#       6609                    6615
#       1013520 t       c       1589755
#       1013527 c       g       1589762
#       1013546 g       -       1589781
#       1288254 -       g       1864488
#       2711293 g       a       3287528
#       2711348 c       t       3287583
#       3390648 -       t       3966883
#       4651064 t       -       5227300
#
############################################################################ 

############################################################################ 
#####                      OPENING FILES
############################################################################ 


### OPEN  THE INDEX FILE (index.tcov) Store it into an array

open(QUAL, $tcov);
	
@coverage = <QUAL>;

close (QUAL);


### OPEN THE UNGAPPED.tcov file (just store it into a filehandle $fh)

open($fh, $tcov2);

### OPEN POSITION LIST (INDEL_info.txt) store into filehandle $fh2

open($fh2, $list);

### OPEN THE INDEL_location file

open(LOCA, "INDEL_location.txt");

@location = <LOCA>;

close (LOCA);


foreach my $loc (@location) {

    @loca = split("\t", $loc);

    if (defined $loca[6]) {

	chomp $loca[6];
    }else {
	
	$loca[6] = " ";
	chomp $loca[5];
    }

    $locat{$loca[0]} = "$loca[5]\t$loca[6]";

}

###########################################################################
##### PARSE THE FILE TO GET THE SNP INFORMATION AND PROCESS
###########################################################################

while ($record = get_next_record($fh2)) {
    
    $refscore = 0;
    $quescore = 0;
    $insscore = 0;
    $covref = 0;
    $qualref = 0;
    $covque = 0;
    $qualque = 0;
    @quality = ();
    @qualityQ = ();
    my $stat = Statistics::Descriptive::Full->new();
    my $statQ = Statistics::Descriptive::Full->new();
    
    ### each record is multiline, split at \n
    ##################################################################
    #####  INDEL-QUE-8:1812926 1 bp ID81 (24538 bp) INS on QUE
    #####  R: 1812906 - 1812946 Q: 2895 - 2936
    #####  R: CCTGGTAGTAATAGGGGAAG-GGCTAGATAAGAAACAATTGG
    #####  Q: CCTGGTAGTAATAGGGGAAGAGGCTAGATAAGAAACAATTGG
    ##################################################################

    my @temp = split("\n", $record);

#    print "$temp[0]\n$temp[1]\n$temp[2]\n$temp[3]\n";


    ### split the first line and second line

    my @temp2 = split(" ", $temp[0]);
    my @temp3 = split(" ", $temp[1]);

    ### grab the position of the SNP on reference and query
    ($indel_id) = ($temp2[0] =~ /^(.*)\:.*$/);
    $indlen = $temp2[1];  # Lenght of INDEL
    $end5R = $temp3[1] + 15;   ## end5 reference
    $end3R = $temp3[3] - 15;   ## end3 reference
    $end5Q = $temp3[5] + 15;   ## end5 query
    $end3Q = $temp3[7] - 15;   ## end3 query
    ($length) = ($temp2[4] =~ /^\((.*)$/);  # length of sequence ID#
 

    print "$temp[0]\n$locat{$indel_id}\n$temp[1]\n$temp[2]\n$temp[3]\n";
    #print "$end5Q\t$end3Q\n";
    
    #print "LENGTH = $length\n";
    
    ###########################################################################
    ### GRAB THE asmbl_id of the QUERY and OPEN THE .tcov FILE
    ###########################################################################
    my $covid; 

    if ($temp2[3] =~ /^ID.*revcom$/) {

	($covid) = ($temp2[3] =~ /^ID(.*)revcom$/);

    }else {

	($covid) = ($temp2[3] =~ /^ID(.*)$/);
    }
    
    
    my $tagname = $tag . "_cov.dir";   ### directory where .tcov files are
    my $tcov3 = "$tagname\/$covid\.tcov";  ## asmbl_id .tcov file

    ### OPENING THE .tcov file

    unless (open(QUAL3, $tcov3)) {
	print "Cannot open $tcov3\n";
	exit(0);
    }


    @coverage3 = <QUAL3>;
    close (QUAL3);
    
    my $t3;
    my @at3;

    ### Eliminate the gap : Create an ungapped .tcov file

    for ($t3 = 0; $t3 < @coverage3; $t3++) {

	if ($coverage3[$t3] =~ /^(\d+)/ && !$at3[$1]) {
	    $at3[$1] = $t3;
	}
    }

    ###########################################################################
    #### TRANFORM THE 6615 REFERENCE COORDINATE TO 6609 COORDINATE
    #### Use subroutine pos_trans
    ###########################################################################
    
    $end5R = pos_trans($end5R);
    $end3R = pos_trans($end3R);

    ###########################################################################
    #### COLLECT THE COVERAGE INFORMATION : TREAT INDEL_REF and INDEL_QUE 
    #### separately. Collect REFERENCE AND QUERY COVERAGE
    ###########################################################################

    if ( $temp2[0] =~ /^INDEL-REF/) {

        #print $temp2[0], "\n";
         
	my $i;
        my $c;

	#### REFERENCE COVERAGE
	print "REF: GBA6615\n";
	for ($i = $end5R-1 ; $i < $end3R ; ++$i) {
	    ### grab the offset from index.tcov
	    my $posoff = $coverage[$i];
	    ### position the cursor to the position offset
	    seek($fh, $posoff, 0);
	    ### Grab the line with subroutine get_line
	    $record2 = get_line($fh);
	    ### format for output (put a * in front of right position)
	    
            my @temp4 = split(" ", $record2);
	    
            my $pos_t = pos_trans2($temp4[0]);
	    
            $record2 = "$pos_t\t$temp4[1]\t$temp4[2]\t$temp4[3]\t$temp4[4]\n";
	    
	    #### Compute the quality score and coverages
	    
	    $lencov = length($temp4[3]);   
	    $covref += $lencov;
	    $qualref += $temp4[2];
	    push(@quality, $temp4[2]);

	    
	    #### Formatting add a * at the position of the INDEL
    
	    if ($i == $end5R + 4) {
		$record2 = "*" . $record2;
	    }else{
		$record2 = " " . $record2;
	    }
            print $record2;
	    
	}
	
	######## CALCULATIONS - STATISTICS
       
	$refcov = $covref/($end3R-$end5R+1);
	$stat->add_data(@quality); 
	$median = $stat->median();
	$refscore = $median/$refcov;
	printf " QUALREF = %.2f  COVREF = %.2f \n ", $refscore, $refcov;

	#### QUERY COVERAGE

	   #### COVERAGE WHERE QUERY IS REVERSE-COMPLEMENT
	print "QUE: $temp2[3]\n";
        if ($temp2[3] =~ /^ID.*revcom$/) {
	    
	    ##### CALCULATE THE EXACT POSITION WITH LENGTH - POSITION....
	    $end5Q = ($length - $end5Q + 1);
	    $end3Q = ($length - $end3Q + 1);
	    	    	    
	    ### go through the coverage in reverse order

	    for ($c = $end5Q+1; $c > $end3Q ; --$c) {
		
		my $tline = $coverage3[$at3[$c-1]];

		## complement bases

		$tline =~ tr/AGTC/TCAG/;
		
		my @temp6 = split(" ", $tline);
		
		my $tline2 = "$temp6[0]\t$temp6[1]\t$temp6[2]\t$temp6[3]\t$temp6[4]\n";
		
		#### Compute the quality score and coverages

		$lenque = length($temp6[3]);
		$covque += $lenque;
	        $qualque += $temp6[2];
		push(@qualityQ, $temp6[2]);

		### Add gap (as many as gap length ($indlen)
		### Add * in front of position before the gap

		if ($c == $end5Q-4) {

		    my $gapline4 = " $c\t-\t--\t----\t---------\n" x  $indlen;
		    print $gapline4;
		    print "*$tline2";

		}else {

		    print " $tline2";
		} 

	    }

	    ######## CALCULATIONS - STATISTICS
	    #$quescore = $qualque/$covque;
	    $quecov = $covque/11;
	    $statQ->add_data(@qualityQ); 
	    $medianQ = $statQ->median();
	    $quescore = $medianQ/$quecov;
	    printf " QUALQUE = %.2f COVQUE = %.2f \n",$quescore, $quecov;

	   #### COVERAGE WHERE QUERY IS IN SAME ORIENTATION AS REFERENCE    

        }else {

	    for ($c = $end5Q-1; $c < $end3Q; ++$c) {

		my $tline = $coverage3[$at3[$c+1]];
		my @temp6 = split(" ", $tline);
		
		my $tline2 = "$temp6[0]\t$temp6[1]\t$temp6[2]\t$temp6[3]\t$temp6[4]\n";
		
		#### Compute the quality score and coverages

		$lenque = length($temp6[3]);
		$covque += $lenque;
	        $qualque += $temp6[2];
		push (@qualityQ, $temp6[2]);

		### Add gap (as many as gap length ($indlen)
		### Add * in front of position before the gap
		
		if ($c == $end5Q + 4) {
		    my $gapline3 = " $c\t-\t--\t----\t---------\n" x  $indlen;
		    print $gapline3;
		    
		    print "*$tline2";
		    
		}else {
		    
		    print " $tline2";
		}
	    }

	    ######## CALCULATIONS - STATISTICS
	    #$quescore = $qualque/$covque;
	    $quecov = $covque/11;
	    $statQ->add_data(@qualityQ); 
	    $medianQ = $statQ->median();
	    $quescore = $medianQ/$quecov;
	    printf " QUALQUE = %.2f COVQUE = %.2f \n", $quescore, $quecov;
	}
	
	### Print record separator
	
	print "//\n";


    ##### INDEL-QUE
	
    }elsif ($temp2[0] =~ /^INDEL-QUE/) {
	
	my $i;
        my $c;
	
	#### REFERENCE COVERAGE
	print "REF: GBA6615\n";
	for ($i = $end5R-1 ; $i < $end3R ; ++$i) {
	    
	    ### grab the offset from index.tcov
	    my $posoff = $coverage[$i];
	    
	    ### position the cursor to the position offset
	    seek($fh, $posoff, 0);
	    
	    ### Grab the line with subroutine get_line
	    
	    $record3 = get_line($fh);
	    
            my @temp5 = split(" ", $record3);
	    
            my $pos_t2 = pos_trans2($temp5[0]);
	    
            $record3 = "$pos_t2\t$temp5[1]\t$temp5[2]\t$temp5[3]\t$temp5[4]\n";
	     		
	    #### Compute the quality score and coverages

	    $lencov = length($temp5[3]);
	    $covref += $lencov;
	    $qualref += $temp5[2];
	    push(@quality, $temp5[2]);

	    ### format for output (put a * in front of right position)
	    ### Add gap

	    if ($i == $end5R + 3) {
		
		$record3 = "*" . $record3;
		print $record3;
                #push (@cov3, $record3);
				
		my $gapline2 = " $pos_t2\t-\t--\t----\t---------\n" x  $indlen;
		print $gapline2;
		# push(@cov3, $gapline);
				
	    }else{
	
		$record3 = " " . $record3;

		print $record3;
	    }
	}	
	
	######## CALCULATIONS - STATISTICS
	#$refscore = $qualref/$covref;
	$refcov = $covref/11;
	$stat->add_data(@quality); 
	$median = $stat->median();
	$refscore = $median/$refcov;
	printf " QUALREF = %.2f  COVREF = %.2f \n ", $refscore, $refcov;
	
	##### QUERY COVERAGE

	#### COVERAGE WHERE QUERY IS REVERSE-COMPLEMENT
	print "QUE: $temp2[3]\n";
	if ($temp2[3] =~ /^ID.*revcom$/) {
	    
	    ##### CALCULATE THE EXACT POSITION WITH LENGTH - POSITION....

	    $end5Q = ($length - $end5Q + 1);
	    $end3Q = ($length - $end3Q + 1);
	    
	    #### Go through coverage in reverse order

	    for ($c = $end5Q+1; $c > $end3Q ; --$c) {
		
		my $tline = $coverage3[$at3[$c-1]];
		
		#### Make complement of coverage		
		$tline =~ tr/AGTC/TCAG/;
		my @temp6 = split(" ", $tline);
		
		my $tline2 = "$temp6[0]\t$temp6[1]\t$temp6[2]\t$temp6[3]\t$temp6[4]\n";
		
		#### Compute the quality score and coverages

		$lenque = length($temp6[3]);
		$covque += $lenque;
	        $qualque += $temp6[2];
		push(@qualityQ, $temp6[2]);

		if ($c == $end5Q-4) {
                    
                    print "*$tline2";
		    
		}else {
		    
		    print " $tline2";
		} 
		
	    }

	    ######## CALCULATIONS - STATISTICS
	    #$quescore = $qualque/$covque;
	    $quecov = $covque/($end5Q-$end3Q+1);
	    $statQ->add_data(@qualityQ); 
	    $medianQ = $statQ->median();
	    $quescore = $medianQ/$quecov;
	    printf " QUALQUE = %.2f COVQUE = %.2f \n", $quescore, $quecov;

	#### COVERAGE WHERE QUERY IS IN SAME ORIENTATION AS REFERENCE 

	}else {

	    for ($c = $end5Q-1; $c < $end3Q; ++$c) {

		my $tline = $coverage3[$at3[$c+1]];
		my @temp6 = split(" ", $tline);
		
		my $tline2 = "$temp6[0]\t$temp6[1]\t$temp6[2]\t$temp6[3]\t$temp6[4]\n";
		
		#### Compute the quality score and coverages

		$lenque = length($temp6[3]);
		$covque += $lenque;
	        $qualque += $temp6[2];
		push(@qualityQ, $temp6[2]);

		if ($c == $end5Q+4) {
                    
                    print "*$tline2";
		    
		}else {

		    print " $tline2";
		} 
		
	    }

	    ######## CALCULATIONS - STATISTICS
	    #$quescore = $qualque/$covque;
	    $quecov = $covque/($end3Q-$end5Q+1);
	    $statQ->add_data(@qualityQ); 
	    $medianQ = $statQ->median();
	    $quescore = $medianQ/$quecov;
	    printf " QUALQUE = %.2f COVQUE = %.2f \n", $quescore, $quecov;
	}
	   
	print "//\n";

    } ## elsif INDEL-QUE

} # while loop

#foreach my $sub (@quality) { 
#    print $sub, "\n";
#}

#close (LIST);

close ($fh);
close ($fh2);

exit;


##############################################################################
###               SUBROUTINES
##############################################################################


sub get_line {


    my ($fh) = @_;

    my $offset;
    my $record = '';
    my ($save_input_separator) = $/;

    $/ = "\n";

    $record = <$fh>;
    $/ = $save_input_separator;

    return $record;

}

sub get_next_record {

    my ($fh) = @_;

    my $offset;
    my $record = '';
    my ($save_input_separator) = $/;

    $/ = "//\n";

    $record = <$fh>;
    $/ = $save_input_separator;

    return $record;

}


sub pos_trans {

    my($pos) = @_;
#   print "SUB POS1: $pos\n";
    my $pos2;

    if ($pos < 1013547) {
	$pos2 = $pos + 576235;
    }elsif ( $pos > 1013546 && $pos < 1288255) {
	$pos2 = $pos + 576234;
    }elsif ( $pos > 1288254 && $pos < 3390649) {
	$pos2 = $pos + 576235;
    }elsif ( $pos > 3390648 && $pos < 4651059) {
	$pos2 = $pos + 576236;
    }elsif ( $pos > 4651058) {
	$pos2 = $pos - 4651058;
    }
#    print "SUB POS1: $pos2\n";
    return $pos2;
    
}


sub pos_trans2 {

 my $pos = shift;
   # print "SUB POS1: $pos\n";
    my $pos2;

    if ($pos > 576235 && $pos < 1589782) {
	$pos2 = $pos - 576235;
    }elsif ( $pos > 1589781 && $pos < 1864489) {
	$pos2 = $pos - 576234;
    }elsif ( $pos > 1864488 && $pos < 3966884) {
	$pos2 = $pos - 576235;
    }elsif ( $pos > 3966883 && $pos < 5227295) {
	$pos2 = $pos - 576236;
    }elsif ( $pos > 0 && $pos < 576236 ) {
	$pos2 = $pos + 4651058;
    }
    #print "SUB POS1: $pos2\n";
    return $pos2;

}



