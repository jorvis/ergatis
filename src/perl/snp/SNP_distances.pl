#! /usr/local/bin/perl -w

use strict;
use lib "/home/jravel/lib";
use BeginPerlBioinfo;
use Getopt::Long;
$| = 1;
 
my $filename;
my $help;
my $HELPTEXT = "

USAGE multifasta_to_fasta.pl -F input filename 

Input filename : SNP_Header.txt file

Input the file as it comes out of the SNP pipeline

HELP: -h (This text)


For questions and report bug please Email : jravel\@tigr.org

";


Getopt::Long::config("no_ignore_case");

my $result = GetOptions ("h|help"     => \$help,
			 "F=s"        => \$filename);

if ($result == 0) {
    die ("Command line parsing failed\n")
    }

if ($help) {

    print $HELPTEXT, "\n";
    exit;
}

if (! defined $filename) {
    print STDERR "You must enter a valid input file\n";
    exit(0);
}

### OPEN TE FILE

my @filedata = get_file_data($filename);

### MAKE ARRAY

my $id; ## The ASMBL_ID from $>ID####revcom-1
my $i = 0; ## Counter
my $distline;  ## New line with distance
my @dist; ## array with new line for output
my $idp;  ## sam as $id but for the previous line
my $dist;  ## Calculated distance between 2 SNPs
my $outputfile = "SNP_distance.txt";

### PARSE THROUGH THE SNP_HEADER.TXT FILE


foreach my $line (@filedata) {

    chomp $line;

    my @temp = split(" ", $line);

    chomp $temp[0];    ### THIS IS >ID###revcom-#
    chomp $temp[4];    ### THIS IS THE SNP POSITION IN THE ASSEMBLY

    if ($i == 0) {     ### FIRST LINE NO CALCULATION TO DO PRINT AS IS with -

	$distline = $line . "\t" . "-";
	push(@dist, $distline);
	++$i;
	next;

    }elsif ($i > 0) {   ### ALL OTHER LINES

	### Extract asmbl_id from >IDasmbl_id-revcom-#
	
	if ($temp[0] =~ /revcom/) {   ### TWO TYPES ONE WITH REVCOM 

	    ($id) = ($temp[0] =~ /^>ID(.*)revcom/);
	
	}else{    ### AND ONE WITHOUT
	    
	    ($id) = ($temp[0] =~ /^>ID(.*)-/);

	}
    
       ### Split the previous line at blank spaces

	my @temp2 = split(" ", $filedata[$i-1]);

	chomp $temp2[0];      
	chomp $temp2[4];
	
	if ($temp2[0] =~ /revcom/) {

	    ($idp) = ($temp2[0] =~ /^>ID(.*)revcom/);
	
	}else{
	    
	    ($idp) = ($temp2[0] =~ /^>ID(.*)-/);

	}

	
	### TEST TO SEE IF SNP BELONG TO THE SAME ASSEMBLY

	if ($id eq $idp) {  ### IF IT DOES CALCULATE THE DISTANCE BETWEENT 2 SNPs
	
	    $dist = $temp[4] - $temp2[4];
	    my $distline = $line . "\t" . $dist;
	    push(@dist, $distline);
	    ++$i;
	    next;
	}elsif ( $id != $idp) {  ### IF IT DOESN'T => NEW ASSEMBLY 
	    
	    ## print new line and then the line with -
	    
	    my $distline = "\n" . $line . "\t" . "-";
	    push (@dist, $distline);
	    ++$i;
	    next;
	}
    
    }
}


### PRINT OUTPUT TO FILE

unless(open(DIST, ">$outputfile")) {
    print "Cannot open $outputfile\n";
}

foreach my $line2 (@dist) {

    print DIST $line2, "\n";
}

close (DIST);


exit;


