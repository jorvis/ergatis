#! /usr/local/bin/perl -w

use strict;
use lib "/home/jravel/lib";
use sub;
use Getopt::Long;
$| = 1;
 
my $con;  ## The .1con REFERENCE
my $filename; ## List of all SNPs in SNP_Header format
my $genome; ## Fasta of the REFERENCE (GBA, pXO1, pXO2...)

my $newline1; 


my $start; ## First 5 bp of the oligo
my $end;  ## Last 5 bp of the oligo
my @keys;   ## Array of all unique SNPs
my @affy;
my $newline2;
my $affyout = "affy.txt";
my @snparray; ## Array of SNP at least 12 bp apart pn either side
my @snpelim;  ## Array of SNP NOT 12 bp apart
my $apart = 5;
my $asmbl_id;
my $coords;
my $newline;
my @window;
my $windowout;
my $windowsize;
my $genomesize;
my %id;  ## Hash of all unique SNPs

my $help;
my $HELPTEXT = "

snpwindow will output a tab delimited of the SNPs distributin along a genome.

USAGE affy.pl -F SNP file -O output file -W windown size -S genome size

-F input file SNP_Header.txt file or a concatenation of all the 
SNP_Header files. The format is identical to the output of the SNP
pipeline.


HELP: -h (This text)


For questions and report bug please Email : jravel\@tigr.org

";


Getopt::Long::config("no_ignore_case");

my $result = GetOptions ("h|help"     => \$help,
			 "O=s"        => \$windowout,
			 "W=s"        => \$windowsize,
			 "S=s"        => \$genomesize,
			 "F=s"        => \$filename);

if ($result == 0) {
    die ("Command line parsing failed\n");
    }

if ($help) {

    print $HELPTEXT, "\n";
    exit;
}

if (! defined $filename) {
    print STDERR "You must enter a valid input file (-F)\n";
    exit(0);
}

#if (! defined $coords) {
 #   print STDERR "You must enter a valid coords file (-C)\n";
 #   exit(0);
#}



### OPEN TE FILE

my @filedata = get_file_data($filename);


### PARSE THE FILE : MAKE HASH TABLE

my $i = 1;
 
foreach my $line (@filedata) {  ### Go through each line of the SNPs file


    chomp $line;

    ## SPLIT EACH FIELD
    
    my @field = split (" ", $line);

    my $refpos = $field[1];   ## The position on the REF is $field[2]
    
    ## CREATE THE HASH TABLE

    if (defined $id{$refpos}) { ## Query the Hash to see if already exist
	                        ## If exist go to next
	next;
    
    }else {
  
	$id{$refpos} = $i;     ## If doesn't exist put it in the hash

	$i++;                  ## Counter goes up
	#print $id{$refpos}, "\n";
	next;
    } 
}


   ## ORDER THE ARRAY OF KEYS

@keys = sort { $a <=> $b } keys %id;


print $keys[2];


my $count = @keys;

print $count;
my $z=0;
my $ct;

while ($z < $genomesize) {
    $ct = 0;
    my $top = $z + $windowsize;

    foreach my $val (@keys) {

	if ($val < $top and $val > $z) {
	    $ct++;
	    next;

	}else{
	    next;
	}
   }

    $newline = $z . "-" . $top . "\t" . $ct;
    push (@window, $newline);

    $z +=  $windowsize;


}


#my $windowout = "snpstats.txt";

unless (open(WINDOW, ">$windowout")) {
    print STDERR "Can't open file!!\n\n";
    exit;
}

foreach my $line2 (@window) { 
    print WINDOW $line2, "\n";
}


exit;
