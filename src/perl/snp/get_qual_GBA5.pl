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
my $tcov2 = "/usr/local/projects/anthrax/jravel/GBA_cov.dir/ungapped.tcov";
my $tcov = "/usr/local/projects/anthrax/jravel/index.tcov";
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
#my $pos
Getopt::Long::config("no_ignore_case");

my $result = GetOptions ("h|help"     => \$help,
			 "P=s"        => \$pos);

my $HELPTEXT = "

USAGE -P position on the genome (6615)


The offset (index.tcov file )is used in this script to retrieve information from the .tcov file.

Perl function tell is used to grab the offset of each line. The seek function in index_get_qual.pl
is used to position the pointer to a particular line for retrieval in a large file.


HELP: -h (This text)
For questions and report bug please Email : jravel\@tigr.org

";




if ($result == 0) {
    die ("Command line parsing failed\n")
    }


if (! defined $pos) {
    print STDERR "You must specify a position (-P)\n";
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




#### DO THE TRANSFORMATION ACCORDING TO SUB_TO_FINAL!!!

    if ($pos < 1013547) {
	$pos = $pos + 576235;
    }elsif ( $pos > 1013546 && $pos < 1288255) {
	$pos = $pos + 576234;
    }elsif ( $pos > 1288254 && $pos < 3390649) {
	$pos = $pos + 576235;
    }elsif ( $pos > 3390648 && $pos < 4651059) {
	$pos = $pos + 576236;
    }elsif ( $pos > 4651058) {
	$pos = $pos - 4651058;
    }


my @cov3;

     
#### CREATE AN ARRAY WITH ONLY 5 LINE BEFORE AND AFTER


    my $i;
    push(@cov3, "REFERENCE: GBA6615\n");
    for ($i = $pos-11 ; $i < $pos+10 ; ++$i) {

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
    print "@cov3";
    


#foreach my $sub (@quality) { 
#    print $sub, "\n";
#}

#close (LIST);

#$end =time;
#my $time3 = $end - $start;
#
#printf "Grabing Info done. Elapsed time %.2fs $time3 \n\n";


#my $time4 = $end - $gstart;

#printf "Process done. Elapsed time %.2fs $time4 \n\n";


#close ($fh);
#close ($fh2);

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



