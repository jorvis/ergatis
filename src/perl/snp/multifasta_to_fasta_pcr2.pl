#! /usr/local/bin/perl -w


use strict;
use TIGR::FASTAiterator;
use TIGR::Foundation;
use TIGR::FASTArecord;
use lib "/home/jravel/lib/";
use BeginPerlBioinfo;
use Getopt::Long;

$| = 1;

my $filename;
my @fastafile;
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



my $newline;

my $tf = new TIGR::Foundation;

if (!defined $tf){
    die ("Bad foundation\n");
}

my @errors;

my $fr = new TIGR::FASTAiterator($tf, \@errors, $filename);

if (!defined $fr){
    die ("Bad reader\n");
}

while ($fr->hasNext){





    my $rec = $fr->next();
    
    my $id = $rec->getIdentifier();
    my $head = $rec->getHeader();
    my $body = $rec->getData();


    if ($id =~ /^BA/) {
	next;
    }

    $id =~ s/^\>//;
    $id =~ s/_\D{3}_\d_\d//;
    my $newline = ">$id\n". print_sequence2($body, 60);
    $id =~ s/^\>//;
    $id =~ s/_\D{3}_\d_\d//;
   
    push(@fastafile, $newline);
    print "THIS IS NEWLINE:  $newline\n";
    

    unless(open(FASTA, ">$id")) {
	print "Cannot open $id!\n\n";
	exit(0);
    }
    
    foreach my $fasta (@fastafile) {
	chomp $fasta;
	print FASTA  $fasta;
    }
    
    @fastafile = ();
    close (FASTA);
    
}


exit;



