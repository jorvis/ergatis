#! /usr/local/bin/perl -w

use strict;
use lib "/home/jravel/lib/";
use BeginPerlBioinfo;
use Getopt::Long;
$| = 1;

### DECLARE VARIABLES

my @files;
my $filename;
my $con;
my $help;
my $HELPTEXT = "

USAGE -F input filename tag -C filename.1con

Input filename Tag : Tag of each folder name

.1con : name of the .1con file in the src folder

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
    print STDERR "You must enter a valid 1con file\n";
    exit(0);
}


my $directory = $filename . "_asmbl_seq.dir";

unless(opendir(DIRECTORY, "$directory")) {
	print "Cannot open directory";
	exit;
	}
	
### READ ALL FILES BUT THE . and .. FILES	
	
@files = grep (!/^\.\.?$/, readdir(DIRECTORY));

closedir(DIRECTORY);

### GO THROUGH EACH FILENAME AND RUN THE SCRIPT

foreach my $file (@files) {	

	if ($file =~ /^ID/) {
	
		chomp $file;
		print $file, "\n\n";

#unless (open READ, $file) {
#print "Cannot open $file\n\n";
#exit;
#}

	#### PIPE THE FILE INTO MUMer
	
	system("~/src/SNP/run-mummer3 ~/src/SNP/$con $directory/$file $file ");
	
	
		my $directory2 = $filename . "_mum_align.dir";

	system("mv $file.* $directory2/.");
		
	
}
}
exit;

