#! /usr/local/bin/perl -w


##################################################################
## This Script takes all the assemblies from 
## the tiling constructed by nucMer and makes 
## the revcom of the one that align to GBA (reference)
## All files are copied or created into a folder call
## GAK_asmbl_seq2.dir
##
##
##
## USAGE: revcom_tiling.pl -F infile

## infile is the file from nucmer which as been parsed by tile.pl
## 
##
##
## These files (in fasta format) are used by Daddy.pl to run the 
## SNP pipeline.
##################################################################

use strict;
use lib "/home/jravel/lib/";
use BeginPerlBioinfo;
use Getopt::Long;

$| =1;

### DECLARE VARIABLES

my $deltafile = "tile.doc";
my $filename;
my @keeper;;
my %listS;
my %listE;
my @temp1;
my @keeper2;

my $help = " USAGE -F filename";



Getopt::Long::config("no_ignore_case");

my $result = GetOptions ("h|help"     => \$help,
			 "F=s"        => \$filename);

if ($result == 0) {
    die ("Command line parsing failed\n")
    }

if (! defined $filename) {
    print STDERR "You must enter a valid folder tag\n";
    exit(0);
}


my $directory = $filename . "_asmbl_seq.dir";

#### READ THE FILE INTO A ARRAY

@keeper = get_file_data($deltafile);

my $i = 0;
my $y = 0;

my %dup = ();

foreach my $set (@keeper) {

    @temp1 = split("\t", $set);

    if ($i < 1) {
	++$i;
	next;

    }else{
	$listS{$i} = $temp1[0];
	$listE{$i} = $temp1[1];
    }
    chomp $temp1[7];
    
    print "$temp1[7] - $temp1[6]\n";
    
    if ($i > 1 and defined $dup{$temp1[7]}) {
	
	print "$temp1[7] is duplicated\n";
	
	++$y;
	
	next;
    }
    
    
    $dup{$temp1[7]} = $temp1[0];
    
    
    if ($temp1[6] =~ /\-/) {
	
	$temp1[7] =~ s/_\D{3}_\d_\d//;
	revcomp($temp1[7]);

    }
    
++$i;

}

print "\n$y\n\n";


exit;

##############################################################
#### SUBROUTINE REVCOMP
##############################################################


sub revcomp {

my ($file) = @_;


### READ IN THE CONTENTS OF THE FILE

my $inputfile =  $file;



my @filedata = get_file_data($inputfile);


### EXTRACT SEQUENCE INFO

my $dna = extract_sequence_from_fasta_data(@filedata);

### CALCULATE THE REVERSE OMPLEMENT OF $DNA

my $revcom = revcom($dna);

### PRINT TO FILE 

my $outputfile = $file . "revcom";

unless (open(FILEOUT, ">$outputfile")) {
    print STDERR "Can't open file!!\n\n";
    exit;
}
print FILEOUT ">reverse-complement of $file\n";
print FILEOUT $revcom;

#foreach my $line2 (@newarray) {
#    print FILEOUT $line2, "\n";
#}

close (FILEOUT);

system("mv $file reverse.dir/$file");

}



