#! /usr/local/bin/perl -w


use strict;
use lib "/home/jravel/lib/";
use BeginPerlBioinfo;
use Getopt::Long;

$| = 1;

my $filename;
my @filename;
my $outputfile;
my $record;

my $help;
my $HELPTEXT = "

USAGE multifasta_to_fasta.pl -F input filename

Input filename : multifasta file

The ouputfile name of each individual sequence will be the name after the > 


HELP: -h (This text)


For questions and report bug please Email : jravel\@tigr.org

";


Getopt::Long::config("no_ignore_case");

my $result = GetOptions ("h|help"     => \$help,
			 "F=s"        => \$filename);
			# "O=s"        => \$outputfile);

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



### OPEN THE MULTIFASTA FILE

open(FH, $filename);

### PARSE MULTIFASTA SET INPUT SEPARATOR TO ">"

$/ = ">";

while (defined ( $record = <FH>)) {

    my @temp = split("\n", $record);
    
    my $count = @temp;
    
    my $i = 1;
    my @fastafile;
    
    foreach my $line (@temp) {
	
	chomp $line;
	
	if ($i == 1) {
	    my @temp2 = split(" ", $line);
	    
	    print $line, "\n";
	    push(@fastafile, ">" . $line);		
	    
	    my @temp2 = split(" ", $line);
	    chomp $temp2[0];
	    $temp2[0] =~ s/^\>//;
	    print $temp2[0], "\n";
	    $outputfile = $temp2[0];
	    print $temp2[0], "\n";
	}else{
	    
	    push(@fastafile, $line);
	}		
	
	++$i;
    }
    
    if ($outputfile eq "") {
	next;
    }else{
	
	unless(open(FASTA, ">$outputfile")) {
	    print "Cannot open $outputfile!\n\n";
	    exit(0);
	}
	
	
	foreach my $fasta (@fastafile) {
	    chomp $fasta;
	    print FASTA  $fasta, "\n";
	}
	
	close (FASTA);
    }
}

exit;



