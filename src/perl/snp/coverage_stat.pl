#! /usr/local/bin/perl -w

use strict;
use lib "/home/jravel/lib/";
use BeginPerlBioinfo;
use Getopt::Long;
$| = 1;
 
my $filename;
my @filename;
my $highcoverage;
my $i;
my $info;
my $data;
my $outputfile;
my @list;
my @reject;
my $help;
my $HELPTEXT = "

USAGE multifasta_to_fasta.pl -F input filename -C high coverage

Input filename : SNP_Header.txt file
high Coverage, the highest coverage in the file

This program will count the number SNP for each coverage 
The highest coverage is entered by the user

HELP: -h (This text)


For questions and report bug please Email : jravel\@tigr.org

";


Getopt::Long::config("no_ignore_case");

my $result = GetOptions ("h|help"     => \$help,
			 "F=s"        => \$filename,
			 "O=s"        => \$outputfile,
			 "C=s"        => \$highcoverage);

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
if (! defined $outputfile) {
    print STDERR "You must enter a valid output filename\n";
    exit(0);
}

if (! defined $highcoverage) {
    print STDERR "You must enter a valid coverage\n";
    exit(0);
}

$info = "Coverage count for $filename";

    push(@filename, $info);


for ($i = 0; $i < $highcoverage + 1; ++$i) {

    open(READ, "grep 'COV: $i ' $filename | wc |");   

my @temp = <READ>;
close (READ);
foreach my $line (@temp) {
    my @temp2 = split(" ", $line);

    $data = $i."X coverage : $temp2[0]";
    push(@filename, $data);
}
}



foreach my $line2 (@filename) { 
    print  $line2, "\n";
}

my @file = get_file_data($filename);

foreach my $snp (@file) {

    my @temp3 = split(" ", $snp);

    if ($temp3[12] > 2 and $temp3[14] > 40) {

	push(@list, $snp);
    } else {
	push(@reject, $snp);
    }

}

$outputfile = $outputfile . ".quality";

unless (open(QUALITY, ">$outputfile")) {
    print STDERR "Can't open file!!\n\n";
    exit;
}

foreach my $sub (@list) { 
   print QUALITY $sub;
}

#print QUALITY "\nCOVERAGE STATISTICS:\n\n";

#foreach my $line3 (@filename) { 
  ##  print QUALITY $line3, "\n";
#}

close (QUALITY);

exit;
