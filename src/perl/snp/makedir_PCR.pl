#! /usr/local/bin/perl -w

#########################################################
## MAKEDIR_PCR is  a script to prepare for ANALYSIS_PCR
## 
## IT DEALS WITH FLAT FILE ASSEMBLY, RETRIEVE INFORMATION 
## FROM FLAT FILE AND GENERATE .CONTIG, .QUAL, and .tcov/

## It creates directories, fill them up with data:
## TAG_asmbl_seq.dir : Sequence of each assembly in fasta
## TAG_contig.dir : Contig file and qual file
## TAG_mum_align.dir : Mummer output
## TAG_cov.dir : Data on Coverage, necessary for GetTotalCoverage
##
## This program takes as input the output of nucmer ran through show-tiling
## This output is taken as is and set as tile.doc in the folder to run the 
## analysis.
##
## First run parse_fasta.pl to eliminate the fasta that are not one assembly
## Once ran, the ANALYSIS_PCR program need to be ran.
## last update from version makedir4 10/29/02
## Added the integration of nucmer and show-tiling
#########################################################



use strict;

use lib "/home/jravel/lib/";
use BeginPerlBioinfo;
use Getopt::Long;
$| = 1;


my $filename; ## FOLDER TAG
my $start; ## COLUMN 1 START POSITION ON REF
my $end; ## COLUMN 2 END POSITION ON REF
my $con; ## .1con file
my $database; ## DATABASE TO BE USED
my @list; 
my $tile = "tile.doc";
my $query; ## multifasta file of all assemblies
my $directory;
### HELP FILE INSTRUCTIONS

my $help;
my $HELPTEXT = "

USAGE -F folder tag -D Directory -Q multifasta -C filename .1con 

Folder tag : All folder create will be tag_.....dir

multifasta: list of assemblies in multifasta

filename .1con: fasta file of the reference in ~/src/SNP/

directory: path to the directory where the assemblyies are for the PCR products

HELP: -h (This text)
For questions and report bug please Email : jravel\@tigr.org

";


### OPTION TO RUN THE PROGRAM

Getopt::Long::config("no_ignore_case");

my $result = GetOptions ("h|help"     => \$help,
			 "F=s"        => \$filename,
			 "D=s"        => \$directory,
			 "Q=s"        => \$query,
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
    print STDERR "You must specify a filename name for the reference\n";
    exit(0);
}
if (! defined $directory) {
    print STDERR "You must specify a directory name for the contig file\n";
    exit(0);
}

### RUN NUCMER WITH THE REFERENCE AND THE MULTIFASTA ASSEMBLIES

system("nucmer ~/src/SNP/$con $query");

system("show-tiling -c out.delta > tile.doc");

system("rm out.*");


### CREATE THE DIRECTORY

   ## TAG_asmbl_seq.dir

my $dir1 = $filename . "_asmbl_seq.dir";

unless (-e $dir1) {
mkdir ($dir1) || die "Cannot mkdir $dir1: $!";
}

   ## TAG_contig.dir

my $dir2 = $filename . "_contig.dir";

unless (-e $dir2) {
mkdir ($dir2) || die "Cannot mkdir $dir2: $!";
}

   ## TAG_mum_align.fir

my $dir3 = $filename . "_mum_align.dir";

unless (-e $dir3) {
mkdir ($dir3) || die "Cannot mkdir $dir3: $!";
}

   ## TAG_cov.dir

my $dir4 = $filename . "_cov.dir";

unless (-e $dir4) {
mkdir ($dir4) || die "Cannot mkdir $dir4: $!";
}



## PARSE THROUGH THE PROCESSED NUCMER OUTPUT

my @tiletolist = get_file_data($tile);
my $i = 0;
my $list = "list.doc"; ## LIST OF CONTIG WITHOUT "ID"


foreach my $mum (@tiletolist) {

    if ($i > 0) { ## SKIP THE FIRST LINE 

	my @mum= split (" ", $mum);

	chomp $mum[7]; ## ID####

	$mum[7] =~ s/_\D{3}_\d_\d//;  ## REMOVE _BAI_1_1 FOR LIST

	push(@list, $mum[7]); ## PUSH INTO ARRAY FOR PRINTING

    }else {

	++$i; ## FIRST THE FIRST LINE

	next;
    }
    ++$i;  ## COUNTER GOES UP
}

### PRINT THE LIST OF ASSEMBLY TO FILE LIST.DOC

unless(open(LIST, ">$list")) {
    print STDERR "Cannot open file $list: $!";
    exit(0);
}

foreach my $idlist (@list) {
    print LIST $idlist, "\n";
}

close (LIST);


### CORE OF THE SCRIPT #########

   ### MOVE list.doc and tile.doc in TAG_asmbl_seq.dir

rename ("list.doc", "$dir1/list.doc") || die "Cannot rename/move the file list.doc: $!";
rename ("tile.doc", "$dir1/tile.doc") || die "Cannot rename/move the file tile.doc: $!";

   ### CHANGE DIRECTORY, NOW IN TAG_asmbl_seq.dir

chdir ($dir1) || die "Cannot cd to $dir1: $!";

#### PROCESS THE list.doc AND DOWNLOAD THE SEQUENCE INTO FASTA FILES

############## USE MULTIFASTA TO FASTA ##############

print "Fetching sequence in $dir1 ...\n\n";

system ("~/src/SNP/multifasta_to_fasta_pcr2.pl -F ../$query");
 
## create individual file with the GBA# as filename and fasta header

  ### CREATE A FOLDER "REVERSE.DIR" FOR ORIGINAL OF REVERSE CONTIGS

unless (-e "reverse.dir") {
mkdir ("reverse.dir") || die "Cannot mkdir reverse.dir: $!";
}
print "Directory reverse.dir created.\n\n";

  ### MOVE list.doc TO TAG_contig.dir

rename ("list.doc", "../$dir2/list.doc") || die "Cannot rename/move the file list.doc: $!";
 
print "File list.doc move into $dir2.\n\n";

  ### CHECK IF CONTIG IN REVERSE OR FORWARD ORIENTATION
    ## IF REVERSE MAKE REVERSE COMPLEMENT AND STORE ORIGINAL 
    ## IN REVERSE.DIR

print "Processing revcom_tiling.pl ...\n\n";

system ("~/src/SNP/revcom_tiling_pcr.pl -F $filename");
  
  ## CHANGE DIRECTORY GO BACK ONE DIRECTORY

chdir ("../") || die "Cannot cd to ../";

### RUNNING MUMMER AND STORE OUTPUT INTO TAG_mum_align.dir

print "Running Mummer...\n\n";


system ("~/src/SNP/Mummy_pcr.pl -F $filename -C $con");


  ## CHANGE DIR GO TO TAG_contig.dir
 
chdir ($dir2) || die "Cannot cd to $dir2: $!";


###### REPLACE WITH SCRIPT TO MOVE FILE FROM ASSEMBLY FOLDER #####
#### MOVE THE .QUAL and .CONTIG files to the TAG_contig.dir

system ("~/src/SNP/contig_pcr.pl -F $filename -D $directory");


## MOVE list.doc to TAG_cov.dir

rename ("list.doc", "../$dir4/list.doc") || die "Cannot rename/move the file list.doc: $!";

chdir ("../$dir4/") || die "Cannot cd to ../$dir4/";

## PULL DOWN COVERAGE INFORMATION AND STORE INTO TAG_cov.dir
## INVOKE GetTotalCoverage

system ("~/src/SNP/get_tcov2.pl -F $filename");



## GO BACK ONE DIRECTORY

chdir ("../") || die "Cannot cd to ../";

unless(opendir(DIRECTORY, "$dir4")) {   
    print "Cannot open directory $dir4";
    exit;
}

my @dirs = grep (!/^\.\.?$/, readdir(DIRECTORY));  

closedir(DIRECTORY);

foreach my $dirlist (@dirs) {

    chomp $dirlist;

    if ($dirlist =~ /^GBA/) {

	my ($tag1, $pre) = ($dirlist =~ /^(GBA\d*)_\D{3}_\d_\w*(\.\w*)/); 

	my $newname = $tag1 . $pre;
	system ("mv $dir4/$dirlist $dir4/$newname");
    }else {
	next;
    }
}

exit;



