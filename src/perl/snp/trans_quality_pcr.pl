#! /usr/local/bin/perl -w


use strict;
use lib "/home/jravel/lib/";
use sub;
use Getopt::Long;
$| = 1;

my $fh2;
my $tcov2;
my $record;
my @trans;
my $help;

Getopt::Long::config("no_ignore_case");

my $result = GetOptions ("h|help"     => \$help,
			 #"T=s"        => \$tag,
			 "F=s"        => \$tcov2);

my $HELPTEXT = qq~

~;


if ($result == 0) {
    die ("Command line parsing failed\n")
    }

if ($help) {
    print $HELPTEXT, "\n";
    exit;
}


if (! defined $tcov2) {
    print STDERR "You must specify a file name (-F)\n";
    exit(0);
}

#if (! defined $tag) {
#    print STDERR "You must specify a folder name (-T)\n";
#    exit(0);
#}


open($fh2, $tcov2);

my @trans2;
my $i = 0;

while ($record = get_next_record($fh2)) {
    
    ++$i;
    if ($i == 1) {
	next;
    }else {
    push (@trans2, ">GBA");
    my @temp = split("\n", $record);
    my $count = @temp;
 #   print "THIS IS COUNT: $count\n";
    if ($count == 4) {
	@trans = @temp[0..1];
    }else {
	@trans = @temp[0..11];
    }
    foreach my $line (@trans) {

	push (@trans2, "$line\n");
	
    }
    
    push (@trans2, "\/\/\n");
}
}

my $output = "SNP_quality2.txt";
open(OUT, ">$output");

foreach my $rec (@trans2) {
    print OUT $rec;

}



exit;


sub get_next_record {

    my ($fh) = @_;

    my $offset;
    my $record = '';
    my ($save_input_separator) = $/;

    $/ = ">GBA";

    $record = <$fh>;
    $/ = $save_input_separator;

    return $record;

}
