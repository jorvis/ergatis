#! /usr/local/bin/perl -w

use strict;

use lib "/home/jravel/lib/";

use Getopt::Long;
$| = 1;


my $filename;
my $start;
my $end;
my $con;
my $database;

my $help;
my $HELPTEXT = "

USAGE -F folderTag  -C filename.1con

HELP: -h (This text)
For questions and report bug please Email : jravel\@tigr.org

";


Getopt::Long::config("no_ignore_case");

my $result = GetOptions ("h|help"     => \$help,
			 "F=s"        => \$filename,			
			 "C=s"        => \$con);

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
if (! defined $con) {
    print STDERR "You must specify a 1con file name\n";
    exit(0);
}


system ("~/src/SNP/Daddy4.pl -F $filename -C $con");


system ("sort -k13nr SNP_Header.txt > SNP_sort_Header.txt");

system ("~/src/SNP/coverage_stat.pl -F SNP_Header.txt -O SNP_cov -C 20");

system( "~/src/SNP/indel.pl -T $filename -C $con > indel.txt");

system( "~/src/SNP/indel_info.pl -T $filename -C $con");

#system("rm queryins.txt");

#system("rm refins.txt");

exit;
