#! /usr/local/bin/perl -w

#############################################################################
###                           get_index_qual.pl
### 
###
###  This script generate a list of offset for every line of the input file
###  The list is used as an index for the index_get_qual.pl
###
###
#############################################################################



use strict;
use lib "/home/jravel/lib/";
use sub;
use Getopt::Long;
use Data::Dumper;

$| = 1;

my $fh;
my @coverage;
my $help;
my $list;
my %index;
my $record;
my $offset;

Getopt::Long::config("no_ignore_case");

my $result = GetOptions ("h|help"     => \$help,
			 "F=s"        => \$list);


my $HELPTEXT = "

USAGE -F .tcov file

This script will create an offset list for each line in the input file.

The offset is used in the script index_get_qual.pl to retrieve information from the .tcov file.

Perl function tell is used to grab the offset of each line. The seek function in index_get_qual.pl
is used to position the pointer to a particular line for retrieval in a large file.


HELP: -h (This text)
For questions and report bug please Email : jravel\@tigr.org

";


if ($result == 0) {
    die ("Command line parsing failed\n")
    }

if ($help) {

    print $HELPTEXT, "\n";
    exit;
}


if (! defined $list) {
    print STDERR "You must specify a file name (-F)\n";
    exit(0);
}



### OPEN  THE UNGAPPED .tcov file
### put the .tcov file into a filehandle ($fh)


unless(open($fh, $list)) {
    print "Cannot open $list\n";
    exit(0);
}



## GET OFFSET OF FIRST LINE

$offset = tell($fh);

## set record separator to end of line "\n"

$/ = "\n";


### GO THROUGH EACH AND PUSH THE OFFSET INTO AN ARRAY

while (defined ( $record = <$fh>)) {

    push(@coverage, "$offset\n");

    ## GET OFFSET OF NEXT LINE

    $offset = tell($fh);
}



#############################################################################
##      PRINT OUTPUT TO index.tcov
#############################################################################


my $outputfile = "index_pXO2.tcov";

unless (open(FILEOUT, ">$outputfile")) {
    print STDERR "Can't open file!!\n\n";
    exit;
}

foreach my $line2 (@coverage) {
    print FILEOUT $line2;
}

close (FILEOUT);

close ($fh);


exit;

