#! /usr/local/bin/perl -w

use strict;
use lib "/home/jravel/lib";
use BeginPerlBioinfo;
use Getopt::Long;
$| = 1;

my $infile = "list.doc";
my $fol;
my $directory;
my $help;
my $HELPTEXT = "

USAGE -F list.doc

HELP: -h (This text)
For questions and report bug please Email : jravel\@tigr.org

";
##############################################################################
### OPTION TO RUN THE PROGRAM
##############################################################################

Getopt::Long::config("no_ignore_case");

my $result = GetOptions ("h|help"     => \$help,
			 "F=s"        => \$fol,
			 "D=s"        => \$directory);

if ($result == 0) {
    die ("Command line parsing failed\n")
    }

if ($help) {

    print $HELPTEXT, "\n";
    exit;
}

if (! defined $directory) {
    print STDERR "You must enter a valid $directory\n";
    exit(0);
}


##############################################################################
### OPEN THE LIST OF SEQ_NAME
##############################################################################
my $folder = $fol . "_PCR";

my $contigfol = $fol . "_contig.dir";

my @list = get_file_data($infile);

print $directory, "\n\n";

chdir ($directory) || die "Cannot cd to $directory: $!";  ### cd in the directory (i.e. BAI)

system ("pwd");



GRAB:foreach my $dir2 (@list) {

    chomp $dir2;
    print $dir2, "\n";

    unless(opendir(DIRECTORY, "$directory/$dir2")) {   
	print "Cannot open directory TEST $dir2";
	
	chdir ("../");
	exit;
    }
    
    my @dirs = grep (!/^\.\.?$/, readdir(DIRECTORY));  
    
    closedir(DIRECTORY);

    my $count = @dirs;
    
    if ($count == 0 ) {
	next;
    }

    chdir ($dir2) || die "Cannot cd back: $!";
    
    my $dir3 = $dir2 ."_PC_" . $count;  ### FOLDER NAME _PC_#

    unless(opendir(DIRECTORY2, "$dir3")) {   
	  print "Cannot open directory $dir2";
	  exit;
      }
	
     ### READ ALL FILES BUT THE . and .. FILES	
      
    my @dirs2 = grep (!/^\.\.?$/, readdir(DIRECTORY2));  
    my $dir4;
    my $out1;
    my $out2;
    foreach my $direct (@dirs2) {
	


	print "THIS IS THE CONTENT: $direct\n";
	
	if ($direct =~ /^asm_2003/) {
	    
	    print "THIS IS IT!!: $direct\n";
	    
	    ($dir4) = ($direct);
	    chomp $dir4;
	}
	if ( $direct =~ /^GBA\d*_\D{3}_\d\.qual/) {

	    $out1 = $direct;
	    chomp $out1;
	}

    }   
    
    
    closedir(DIRECTORY2);
    
  
    ($out2) = ($out1 =~ /^(GBA\d*_\D{3}_\d)/);    
    
    chdir ($dir3) || die "Cannot move into $dir3: $!";
    
    my ($out) = ($out1 =~ /^(GBA\d*)_\D{3}.*/);

    my $dotqual = $out . ".qual";
    
    system ("cp $out1 /usr/local/projects/$fol/jravel/$folder/$contigfol/$dotqual");

    chdir ("$dir4") || die "Cannot cd back3: $!";

    

    my $input = $out2 . ".fasta";

      
    print "THIS IS INPUT: $input\n";

    #my ($newline, $ctr) = tofasta($input);
      
   # if ($ctr > 1) {
	#chdir ("../../../") || die "Cannot restart to next line!: $!";
	#system ("rm /usr/local/projects/$fol/jravel/$folder/$contigfol/$dotqual");
	#next GRAB;
    #}
    
    my $dir5 = $out2 . ".align";
    
    unless(opendir(DIRECTORY3, "$dir5")) {   
	print "Cannot open directory $dir5";
	  
	chdir ("../../../") || die "Cannot restart to next line!: $!";
	
	next GRAB;
    }
    
    my @dirs3 = grep (!/^\.\.?$/, readdir(DIRECTORY3));  
    
    my $dotcontig = $out . ".contig";
    foreach my $direct (@dirs3) {

	system("cp  $dir5/$direct  /usr/local/projects/$fol/jravel/$folder/$contigfol/$dotcontig");
	
	
    }   
    
    
    closedir(DIRECTORY3);
    
    
    
    chdir ("../../../") || die "Cannot cd back: $!";
    
}

exit;



sub tofasta {

use TIGR::FASTAiterator;
use TIGR::Foundation;
use TIGR::FASTArecord;
use lib "/home/jravel/lib/";
use BeginPerlBioinfo;

my ($asmb) = @_;

print "THIS IS ASMB; $asmb\n";

my $newline;

my $tf = new TIGR::Foundation;

if (!defined $tf){
    die ("Bad foundation\n");
}

my @errors;

my $fr = new TIGR::FASTAiterator($tf, \@errors, $asmb);

if (!defined $fr){
    die ("Bad reader\n");
}
my $line;
my $i = 0;
while ($fr->hasNext){
    my $rec = $fr->next();
    
    my $id = $rec->getIdentifier();
    my $head = $rec->getHeader();
    my $body = $rec->getData();
    my $newline = "$head\n". print_sequence2($body, 60);
    print "THIS IS NEWLINE:  $newline\n";
    $line .= $newline;
    $i++;

}
return ($line, $i);



}
