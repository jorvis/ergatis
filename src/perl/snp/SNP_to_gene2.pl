#! /usr/local/bin/perl -w

#########################################################
#### LIST THE GENE WITH SNP                     #########
#### Given the output SNP_Header.txt list all the gene
#### and the com_name in a tabulated list 
#########################################################


use strict;
use lib "/home/jravel/lib/";
use BeginPerlBioinfo;

use strict;
use Getopt::Long;

$| = 1;


### declare variables

my $line;
my $newline;
my $row;
my @coordarray = ();      ### Array of ORF# END5 END3
my @matcharray = ();      ### Array of MOTIF POSITION
my @SNP_to_gene = ();
my @temp2 = ();
my @temp3 = ();
my @temp = ();            ### tempory array to store line of @coordarray
my $outputfile;
my $row2;
my $filename1;
my $filename2;
my $tag;
my $help;
my $HELPTEXT = "USAGE -F coordinate filename -R SNP list file name\n";


Getopt::Long::config("no_ignore_case");

my $result = GetOptions ("h|help"     => \$help,
			 "F=s"        => \$filename1,
			 "R=s"        => \$filename2,
			 "T=s"        => \$tag);

if ($result == 0) {
    die ("Command line parsing failed\n")
    }
if ($help) {

    print $HELPTEXT, "\n";
    exit;
}

if (! defined $filename1) {
    print STDERR "You must enter a valid input file\n";
    exit(0);
}
if (! defined $filename2) {
    print STDERR "You must specify a reference file name\n";
    exit(0);
}
if (! defined $tag) {
    print STDERR "You must specify a Tag name for output file\n";
    exit(0);
}


$outputfile = $tag . "_to_gene.txt";

#### READ THE REFERENCE FILE INTO AN ARRAY

@coordarray = get_file_data($filename1);

@matcharray = get_file_data($filename2);


foreach $line (@coordarray) {
    chomp $line;
    @temp = split ("::", $line);
    

    my $end5 = $temp[1];
    my $end3 = $temp[2];
    
    if ($end3 > $end5) {
	
	foreach $row (@matcharray) {
	    chomp $row;
	    @temp2 = split(' ', $row);

	    
	    if ($temp2[2] > $end5 and $temp2[2] < $end3) {
		push (@SNP_to_gene, $temp2[0] . "\t" . $temp2[2] ."\t" .$temp2[3]. "\t" . $temp2[4] . "\t" . $temp[0] . "\t" . $temp[3]);
	       
		next;
	    }
	}
    }elsif ($end5 > $end3) {

 	foreach $row2 (@matcharray) {
	    chomp $row2;

	    @temp3 = split(' ', $row2);

	    if ($temp3[2] > $end3 and $temp3[2] < $end5) {
		push (@SNP_to_gene,  $temp3[0] . "\t" . $temp3[2] . "\t" .$temp3[3]. "\t" . $temp3[4] . "\t" .  $temp[0] . "\t" . $temp[3]);
	       
		next;
	    }
	}

    }
}

### PRINT TO FILE 

unless (open(FILEOUT, ">$outputfile")) {
    print STDERR "Can't open file!!\n\n";
    exit;
}

foreach my $line2 (@SNP_to_gene) {
    print FILEOUT $line2, "\n";
}

close (FILEOUT);

exit;





exit;


##########################################################
#### SUBROUTINE PRINT TABLE
##########################################################

sub printTable {
    my $fileout = shift; 
    my $rt = shift;
    my %ht = %$rt;
    my @keys = sort keys %ht;
    
    open(OUT, ">>$fileout") or die "Failed to open $fileout ($!)";

    foreach my $i (@keys){
	print "$i\t" . $ht{$i} . "\n";
	print OUT "$i\t" . $ht{$i} . "\n"  or die ("Failed to write to $fileout ($!)");
    }
	print "\n\n";
    print OUT "\n\n\n";
    close(OUT) or die ("Failed to close $fileout ($!)");
}


    


