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
use Getopt::Long;
$| = 1;

my $start;
my $gstart;
my $end;
my $pos;
my $help;
my $newcov;
my @quality;
my $tcov2 = "/usr/local/projects/anthrax/jravel/pXO2_cov.dir/ungapped_pXO2.tcov";
my $tcov = "/usr/local/projects/anthrax/jravel/index_pXO2.tcov";
my $add;
my $list;
my @coverage;
my $fh;
my @list;
my $fh2;
my $record;
my @location;
my @loca;
my %locat;
my %quali;
my $qualine;
my $fh3;
my $refscore;
my $refqual;
my $refcov;

Getopt::Long::config("no_ignore_case");

my $result = GetOptions ("h|help"     => \$help,
			 "F=s"        => \$list);

my $HELPTEXT = "

USAGE -F list file (SNP_quality.txt)

The list file is the output of the SNP pipeline and looks like this:

>ID1052-1 748158 a-->c 2360
 2355    T       129     TTT     47:47:35
 2356    T       95      TTT     44:22:29
 2357    A       90      AAA     47:25:18
 2358    T       77      TTT     47:13:17
 2359    C       60      CC-     47:13:0
*2360    C       6       ACC     40:10:25
 2361    A       68      AAA     40:9:19
 2362    C       93      CCC     44:15:34
 2363    T       105     TTT     45:28:32
 2364    G       130     GGG     44:49:37
 2365    A       96      AAA     31:40:25

>ID994-1 753325 c-->a 2628
 2623    A       66      AAA     24:24:18
 2624    G       62      GGG     22:18:22
 2625    A       56      AAA     19:17:20
 2626    A       52      AAA     15:17:20
 2627    A       55      AAA     18:17:20
 2627    -       18      C--     18:0:0
*2628    A       43      AAA     18:12:13
 2629    C       12      ACC     18:14:16
 2630    A       67      AAA     21:19:27
 2631    A       65      AAA     21:19:25
 2632    G       52      GGG     15:19:18
 2633    A       66      AAAA            15:18:17:16
 ...

The output is the same file with the cutout of the .tcov for that position added to each record.

The offset (index.tcov file )is used in this script to retrieve information from the .tcov file.

Perl function tell is used to grab the offset of each line. The seek function in index_get_qual.pl
is used to position the pointer to a particular line for retrieval in a large file.


HELP: -h (This text)
For questions and report bug please Email : jravel\@tigr.org

";




if ($result == 0) {
    die ("Command line parsing failed\n")
    }


if (! defined $list) {
    print STDERR "You must specify a file name (-F)\n";
    exit(0);
}


##### DO THE TRANSFORMATION FROM 6615 to 6609 ACCORDING TO:

############################################################################
# asmbl_id    asm_lend    asm_rend    sub_asmbl_id sub_asm_lend sub_asm_rend
# ----------- ----------- ----------- ------------ ------------ ------------
#       6615           1     4651062         6609       576236      5227297
#       6615     4651063     5227297         6609            1       576235
############################################################################ 




### OPEN  THE INDEX FILE (index.tcov) Store it into an array

open(OUT, ">SNP_final.txt");


$start = time;
$gstart = $start;


open(QUAL, $tcov);
	
@coverage = <QUAL>;

close (QUAL);

    ## time the opening of index and storing into an array

$end =time;
my $time1 = $end - $start;

printf "Opening .tcov file done. Elapsed time %.2fs  \n\n", $time1;

$start = time;

### OPEN THE UNGAPPED.tcov file (just store it into a filehandle $fh)

open($fh, $tcov2);

## time the storing in $fh (should be 0)

$end =time;
my $time2 = $end - $start;

printf "Putting ungapped.tcov in filehandle done. Elapsed time %.2fs \n\n", $time2;


### OPEN POSITION LIST (SNP_quality.txt) store into filehandle $fh2

$start = time;

open($fh2, $list);

### time the storing into $fh2

$end =time;
$time2 = $end - $start;

printf "Opening list file done. Elapsed time %.2fs  \n\n", $time2;

### OPEN THE INDEL_location file

open(LOCA, "SNP_location.txt");

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
#### OPEN THE SNP_align2.txt file

open($fh3, "SNP_align2.txt");

while ($qualine = get_next_record($fh3)) {

    my @qual2 = split ("\n", $qualine);

    my @qual3 = split(" ", $qual2[0]);

    my $id2 = $qual3[0];

    my $quescore = $qual3[12]/$qual3[10];
    $quescore = int($quescore*100)/100;

    $quali{$id2} = "$qual3[0] $qual3[1] $qual3[2] $qual3[3] $qual3[4] $qual3[9] $qual3[10] $qual3[11] $qual3[12] QUAL: $quescore\n$qual2[1]\n$qual2[2]\n$qual2[3]";

}



##### PARSE THE LIST TO GET THE SNP POSITION AND PROCESS


while ($record = get_next_record($fh2)) {

    ### each record is multiline, split at \n
    my @temp = split("\n", $record);

    ### split the first line
    my @temp2 = split(" ", $temp[0]);
    
    ### grab the position of the SNP on reference 
    $pos = $temp2[1];
    
    my $rec2 =  "$temp[1]\n$temp[2]\n$temp[3]\n$temp[4]\n$temp[5]\n$temp[6]\n$temp[7]\n$temp[8]\n$temp[9]\n$temp[10]\n$temp[11]\n";

    if ($temp2[0] =~ /^\>ID\d+\-\d+/) {
	
	$rec2 =~ tr/ATGC/TACG/;
	my @temp3 = split("\n", $rec2);
	my $rec3 = "$quali{$temp2[0]}\n$locat{$temp2[0]}\n$temp3[10]\n$temp3[9]\n$temp3[8]\n$temp3[7]\n$temp3[6]\n$temp3[5]\n$temp3[4]\n$temp3[3]\n$temp3[2]\n$temp3[1]\n$temp3[0]";
	print OUT $rec3, "\n";

	push(@quality, $rec3);

    }elsif ($temp2[0] =~ /revcom/) {
	
	$rec2 =  "$quali{$temp2[0]}\n$locat{$temp2[0]}\n$temp[1]\n$temp[2]\n$temp[3]\n$temp[4]\n$temp[5]\n$temp[6]\n$temp[7]\n$temp[8]\n$temp[9]\n$temp[10]\n$temp[11]\n";


	print OUT $rec2;

	push(@quality, $rec2);
	
    }

# DO THE POS TRANSFORMATION


    if ($pos < 51578) {
	$pos = 51578 - $pos ;
	
    }elsif ($pos > 51577) {
	$pos = (94830 - $pos) + 51577 ;
	
    }

    my @cov3;
     
#### CREATE AN ARRAY WITH ONLY 5 LINE BEFORE AND AFTER

    my $i;

    push(@cov3, "REFERENCE: GBX17224\n");

    for ($i = $pos-6 ; $i < $pos+5 ; ++$i) {

	### grab the offset from index.tcov
	my $posoff = $coverage[$i];
     
	### position the cursor to the position offset
	seek($fh, $posoff, 0);
	
	### Grab the line with subroutine get_line

	$record = get_line($fh);

	### format for output (put a * in front of right position)

	my @temp4 = split(" ", $record);

	my $pos_t = pos_trans2($temp4[0]);

	$record = "$pos_t\t$temp4[1]\t$temp4[2]\t$temp4[3]\t$temp4[4]\n";

	if ($i == $pos-1) {

	    $refcov = length($temp4[3]);
	    $refscore = $temp4[2]/$refcov;
	    $refscore = int($refscore*100)/100;
	    $refqual = $temp4[2];
	    $record = "*" . $record;

	}else{
	
	    $record = " " . $record;
	}
	
	push (@cov3, $record);
	
    
    }
    push (@cov3, "COV: $refcov CBQUAL: $refqual QUAL: $refscore\n\/\/\n");
    print OUT "@cov3";
    
}


$end =time;
my $time3 = $end - $start;

printf "Grabing Info done. Elapsed time %.2fs $time3 \n\n";


my $time4 = $end - $gstart;

printf "Process done. Elapsed time %.2fs $time4 \n\n";


close ($fh);
close ($fh2);
close (OUT);
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

    $/ = "\/\/\n";

    $record = <$fh>;
    $/ = $save_input_separator;

    return $record;

}


sub pos_trans2 {

 my $pos = shift;

    my $pos2;
 
 if ($pos > 51577 && $pos  < 94830) {

     $pos2 = (94830 - $pos) + 51577;
	 
 }elsif ($pos < 51577) {
 
    $pos2 = 51578 - $pos ;
 
 }


 return $pos2;

}
 
