#! /usr/local/bin/perl -w

use strict;

use lib "/home/jravel/lib/";
use BeginPerlBioinfo;
use Getopt::Long;
$| = 1;

my $coords; ## _com.coords file
my $tag; ## FOLDER TAG
my $genome;
my $SNPfile1 = "SNP_Header.txt";
my $SNPfile2 = "SNP_sort_Header.txt";
my $INDELfile = "INDEL_Header.txt";

### HELP FILE INSTRUCTIONS

my $help;
my $HELPTEXT = qq~

USAGE -C _com.coords -T Tag (for output) -G G, 1, 2

-C _com.coords : .coords file with com_name added

-T Tag : a Tag for the output.

-G chromosome : "G", pXO1: "1", pXO2: "2"

This script assumes that the makedir and analysis script have been run, and 
SNP_sort_header.txt, SNP_align.txt, SNP_quality.txt and SNP_Header.txt are present in the directory the script 
is ran from.

HELP: -h (This text)
For questions and report bug please Email : jravel\@tigr.org

~;


### OPTION TO RUN THE PROGRAM

Getopt::Long::config("no_ignore_case");

my $result = GetOptions ("h|help"     => \$help,
			 "C=s"        => \$coords,
			 "T=s"        => \$tag,
			 "G=s"        => \$genome);

if ($result == 0) {
    die ("Command line parsing failed\n");
    }

if ($help) {

    print $HELPTEXT, "\n";
    exit;
}

if (! defined $coords) {
    print STDERR "You must enter a valid coords file (-C) \n";
    exit(0);
}
if (! defined $tag) {
    print STDERR "You must enter a Tag name (-T)\n";
    exit(0);
}

if (! defined $genome) {
    print STDERR "You must enter a genome type chromosome: G, pXO1: 1, pXO2: 2 (-G)\n";
    exit(0);
}

if ($genome eq "G" or $genome eq  "1" or  $genome eq "2") {
    
}else {
    print STDERR "You must enter a genome type chromosome: G, pXO1: 1, pXO2: 2 (-G)\n";
    exit(0);
}

system("~/src/SNP/SNP_loc.pl -F $SNPfile1 -C $coords -P 2 -O SNP");

system("~/src/SNP/SNP_loc.pl -F $INDELfile -C $coords -P 1 -O INDEL");

#system("~/src/SNP/SNP_to_gene2.pl -F $coords -R $SNPfile2 -T $tag");

#system("~/src/SNP/SNP_intergenic2.pl -P $SNPfile1 -O $coords -T $tag");

system("~/src/SNP/SNP_distance_pcr.pl -F $SNPfile1");

print "Transforming the file .... align.txt ";

system("~/src/SNP/trans_align_pcr.pl -F SNP_align.txt");

print "...quality.txt\n\n";

system("~/src/SNP/trans_quality_pcr.pl -F SNP_quality.txt");


print "\n Getting coverage information for the reference...  ";

if ($genome eq "G") {

    system ("~/src/SNP/index_get_qual_pcr.pl -F SNP_quality2.txt");
    print "Done.\n\n";
    print "\n Getting information on Synonymous-NonSynonymous SNPs...\n ";
    system ("~/src/SNP/syn_nonsyn_pcr.pl -F SNP_final.txt -T $tag");
    print "\nSynonymous-NonSynonymous SNPs - Done.\n\n";

}elsif ($genome eq "1") {

    system ("~/src/SNP/index_get_qual_pXO1.pl -F SNP_quality2.txt");
    print "Done.\n\n";
    print "\n Getting information on Synonymous-NonSynonymous SNPs...\n ";
    system ("~/src/SNP/syn_nonsyn_pXO1.pl -F SNP_final.txt -T $tag");
    print "\nSynonymous-NonSynonymous SNPs - Done.\n\n";

}elsif ($genome eq "2") {

    system ("~/src/SNP/index_get_qual_pXO1.pl -F SNP_quality2.txt");
    print "Done.\n\n";
    print "\n Getting information on Synonymous-NonSynonymous SNPs...\n ";
    system ("~/src/SNP/syn_nonsyn_pXO2.pl -F SNP_final.txt -T $tag");
    print "\nSynonymous-NonSynonymous SNPs - Done.\n\n";

}

print "Running SNP parser and Gathering summary information...\n";

system ("~/src/SNP/STATS.dir/SNP_parser_pcr.pl -F output.final");
system ("~/src/SNP/STATS.dir/final_parse.pl -F SNPs2.txt");
system ("~/src/SNP/STATS.dir/stats_pcr.pl -F HL.snp -T $tag");


print "\n\nSNP Process Done.\n\n";

print "Now processing INDEL...\n\n";

system ("~/src/SNP/cov_indel4.pl -T $tag -F INDEL_info.txt");

print "INDEL Process Done.\n\n";


exit;

