#! /usr/local/bin/perl -w


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
my $tcov = "/usr/local/projects/anthrax/jravel/GBA_cov.dir/6609.tcov";


Getopt::Long::config("no_ignore_case");

my $result = GetOptions ("h|help"     => \$help,
			 "P=s"        => \$pos);

if ($result == 0) {
    die ("Command line parsing failed\n")
    }


if (! defined $tcov) {
    print STDERR "You must specify a coverage file name\n";
    exit(0);
}


##### DO THE TRANSFORMATION FROM 6615 to 6609 ACCORDING TO:

############################################################################
# asmbl_id    asm_lend    asm_rend    sub_asmbl_id sub_asm_lend sub_asm_rend
# ----------- ----------- ----------- ------------ ------------ ------------
#       6615           1     4651062         6609       576236      5227297
#       6615     4651063     5227297         6609            1       576235
############################################################################ 


my $add;

if ($pos < 4651063) {
    $pos = $pos + 576235;
    print $pos, "\n";
    $add = 1;
}elsif ($pos > 4651062) {
    $pos = $pos - 4651062;
    print $pos, "\n";
    $add = 0;
}


### ELIMINATE THE GAPS AND DUPLICATE POSITION SO LINE# = POSITION NUMBER
$start = time;
$gstart = $start;
my @coverage;

open(QUAL, $tcov);
	
@coverage = <QUAL>;
close (QUAL);

$end =time;
my $time1 = $end - $start;

printf "Opening file done. Elapsed time %.2fs  \n\n", $time1;


$start = time;
my $t;
my @at;
my @cov4;
for ($t = 0; $t < @coverage; $t++) {
    
    if ($coverage[$t] =~ /^(\d+)/ && !$at[$1]) {
	$at[$1] = $t;
	push (@cov4, $coverage[$t]);
    }
}


$end =time;
my $time2 = $end - $start;

printf "Parsing .tcov file done. Elapsed time %.2fs \n\n", $time2;

my @cov2;
my @cov3;



     
#### CREATE AN ARRAY WITH ONLY 10 LINE BEFORE AND AFTER

$start = time;

 @cov3 = @coverage[$at[$pos-10]..$at[$pos+10]];

 foreach my $covline (@cov3) { ### GO THROUGH EACH LINE
		
		my @cov2 = split(" ", $covline);

		if ($cov2[0] == $pos) {  ### TAG THE LINE FOR THE SNP WITH *
		
		    if ($add == 1) {
			$cov2[0] = $cov2[0] - 576235;
		    }elsif ($add == 0) {
			$cov2[0] = $cov2[0] + 4651062;
		    }
		    
		    $newcov = "*" . "$cov2[0] \t $cov2[1] \t $cov2[2] \t $cov2[3]   \t $cov2[4] \t $cov2[5]";
		    push (@quality, $newcov);

		}else {  ### BASICALLY DO NOTHING HERE BUT ADD THE TAB

		     if ($add == 1) {
			$cov2[0] = $cov2[0] - 576235;
		    }elsif ($add == 0) {
			$cov2[0] = $cov2[0] + 4651062;
		    }

		    $newcov = " " . "$cov2[0] \t $cov2[1] \t $cov2[2] \t $cov2[3]   \t $cov2[4] \t $cov2[5]";
		    push (@quality, $newcov);
		}
	    }  

foreach my $sub (@quality) { 
    print $sub, "\n";
}

$end =time;
my $time3 = $end - $start;

printf "Grabing Info done. Elapsed time %.2fs $time3 \n\n";


my $time4 = $end - $gstart;

printf "Process done. Elapsed time %.2fs $time4 \n\n";

### PRINT TO FILE 
my $outputfile = "ungapped.tcov";

unless (open(FILEOUT, ">$outputfile")) {
    print STDERR "Can't open file!!\n\n";
    exit;
}

foreach my $line2 (@cov4) {
    print FILEOUT $line2;
}

close (FILEOUT);



exit;

