#! /usr/local/bin/perl -w


use strict;
use lib "/home/jravel/lib/";
use sub;
use Getopt::Long;
$| = 1;

my $orientation;
my $tcov;

my $length;
my $pos;
my $help;
my $newcov;
my @quality;


Getopt::Long::config("no_ignore_case");

my $result = GetOptions ("h|help"     => \$help,
			 "O=s"        => \$orientation,
			 "L=s"        => \$length,
			 "P=s"        => \$pos,
			 "C=s"        => \$tcov);

if ($result == 0) {
    die ("Command line parsing failed\n")
    }
if (! defined $orientation) {
    print STDERR "You must enter a valid orientation parameter\n";
    exit(0);
}

if (! defined $length) {
    print STDERR "You must enter a valid length\n";
    exit(0);
}

if (! defined $tcov) {
    print STDERR "You must specify a coverage file name\n";
    exit(0);
}


my @coverage;

open(QUAL, $tcov);
	
@coverage = <QUAL>;
close (QUAL);

my $t;
my @at;
my @cov4;

for ($t = 0; $t < @coverage; $t++) {

    if ($coverage[$t] =~ /^(\d+)/ && !$at[$1]) {
	$at[$1] = $t;
	push (@cov4, $coverage[$t]);
    }
}

my @cov2;
my @cov3;
my $posM;

if ($orientation eq "R"){  ## IF THE FASTA FILE CONTAIN revcom THEN COMPUTE COORDINATE
	    $posM = ($length - $pos + 1);  ### LENGTH OF ASSEMBLY - SNP pos + 1 (start at 0)

	}else{  ## IF ON FORWARD STRAND JUST KEEP THE POSITION OF THE SNP
	    $posM = $pos;
        }


 @cov3 = @coverage[$at[$posM-15]..$at[$posM+15]];  ### CHANGE HERE for more coverage

 foreach my $covline (@cov3) { ### GO THROUGH EACH LINE
		
		my @cov2 = split(" ", $covline);

		
		if ($cov2[0] == $posM) {  ### TAG THE LINE FOR THE SNP WITH *
		    
		 #   $pval = $cov2[2];
		#    $lval = length($cov2[3]);
		    
		    $newcov = "*" . "$cov2[0] \t $cov2[1] \t $cov2[2] \t $cov2[3]   \t $cov2[4]";
		    push (@quality, $newcov);

		}else {  ### BASICALLY DO NOTHING HERE BUT ADD THE TAB

		    $newcov = " " . "$cov2[0] \t $cov2[1] \t $cov2[2] \t $cov2[3]   \t $cov2[4]";
		    push (@quality, $newcov);

		}

	    }  


foreach my $sub (@quality) { 
    print $sub, "\n";
}

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

