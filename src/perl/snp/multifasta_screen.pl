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

## GO THROUGH EACH FASTA FILE 

	
	my @fastafile;

while (defined ( $record = <FH>)) {

        ## SPLIT THE RECORD AT \n
	my @temp = split("\n", $record);
	
	my $count = @temp;

	my @temp3 = split(" ", $temp[0]);
	my $i = 1;
	if ($temp3[1] < 1400) {

	    next;
	}else{
		foreach my $line (@temp) {
				
			chomp $line;

			if ($i == 1) {

					
			### NAME OF THE OUTPUT FILE
			    ## SPLIT THE LINE 
			my @temp2 = split(" ", $line);
			chomp $temp2[2];
			    ## DELETE > FROM THE NAME
		#	$temp2[0] =~ s/^\>//;
		#	print $temp2[0], "\n";
			    ## SET THE OUTPUTFILE TO THE NAME
		#	$outputfile = $temp2[0];
		#	print $temp2[0], "\n";
			
		#	       if ($temp2[2] > 1200) {
	  	 
			
			
			### PUSH FIRST LINE IN THE ARRAY WITH > IN FRONT
			my $line2 = ">" . $line . "\n";
				   print $line2, "\n";
				   push(@fastafile, $line2);
			#       }else{
				#   last;
			      # }

			}else{

			push(@fastafile, $line . "\n");

			}		

		++$i;
		}
	    }
    }
	#	if ($outputfile eq "") {
	#	next;
	#	}else{

$outputfile = $filename . "2";

		unless(open(FASTA, ">$outputfile")) {
    		print "Cannot open $outputfile!\n\n";
    		exit(0);
		}


		foreach my $fasta (@fastafile) {
    		chomp $fasta;
    		print FASTA  $fasta;
		}

		close (FASTA);
		


exit;








