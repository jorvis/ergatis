#! /usr/local/bin/perl -w

use strict;
use lib "/home/jravel/lib";
use BeginPerlBioinfo;
use Getopt::Long;
$| = 1;

my $database;
my $help;
my $HELPTEXT = "

USAGE -D database

HELP: -h (This text)
For questions and report bug please Email : jravel\@tigr.org

";


Getopt::Long::config("no_ignore_case");

my $result = GetOptions ("h|help"     => \$help,			
			 "D=s"        => \$database);

if ($result == 0) {
    die ("Command line parsing failed\n")
    }

if ($help) {

    print $HELPTEXT, "\n";
    exit;
}

if (! defined $database) {
    print STDERR "You must enter a valid $database\n";
    exit(0);
}



my $cont;
my $input= "list.doc";
my @temp = get_file_data($input);

foreach $cont (@temp){

    chomp $cont;

    system ("echo $cont > $cont");

    system ("pull_contig -D $database  -A $cont -o $cont");

    system ("rm *.seq");

}

exit;



