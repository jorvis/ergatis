#! /usr/local/bin/perl -w

use strict;
use lib "/home/jravel/lib";
use BeginPerlBioinfo;
use Getopt::Long;
$| = 1;

my $filename;
my $help;
my $HELPTEXT = "

Run the getTotalCoverage Script for the list of contig

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
    print STDERR "You must enter a valid $filename\n";
    exit(0);
}




my $cont;
my $input= "list.doc";
my @temp = get_file_data($input);
my $directory = $filename . "_contig.dir";
foreach $cont (@temp){

    chomp $cont;
   
    print "Getting tcov files for contig $cont...\n\n";

    system ("getCoverage -t -o .  ../$directory/$cont");

    
}

exit;



