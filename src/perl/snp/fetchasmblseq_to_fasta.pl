#! /usr/local/bin/perl -w

############################################################
######### Database interface
############################################################

use strict;
use English;
use Getopt::Long;
use DBI;
use lib "/home/jravel/lib/";
use utilities;
use BeginPerlBioinfo;

$| = 1;


#my($revision) = '$Revision$';
my ($dbserver) = "SYBTIGR";
my ($dbuser) = "access";
my ($dbpassword) = "access";
my ($help);
my ($db);
my ($dbtype) = "Sybase";
my ($version);
my ($fileinput) = "list.doc";
my (@newarray)=();
my ($outputfile);
Getopt::Long::config("no_ignore_case");

my ($result) = GetOptions ("h|help"    => \$help,
			   "D=s"       => \$db,
			   "S=s"       => \$dbserver,
			   "U=s"       => \$dbuser,
			   "P=s"       => \$dbpassword,
                           "F=s"       => \$fileinput);

if ($result == 0) {
    die ("Command line parsing failed\n")
    }

if (! defined $db) {
    print STDERR "You must specify a database\n";
    exit(0);
}
if (! defined $fileinput) {
    print STDERR "You must specify a filename\n";
    exit(0);
}


my ($dbproc) = DBI->connect("dbi:$dbtype:server=$dbserver;packetSize=8092",
			    $dbuser, $dbpassword);

    if (! defined $dbproc) {
	die ("Connection to server $dbserver failed: $!\n");
    }

	$dbproc->do("use $db")  || die ("use $db error: " . $dbproc->errstr);
	$dbproc->do("set textsize 50000000")
	    || die ("set textsize error: " . $dbproc->errstr);


##############################################################
### Getting com_name for each feat_name
##############################################################

my @list = get_file_data ($fileinput);

foreach my $line (@list) {

    chomp $line;

    my @temp = split(' ', $line);


my ($query) = qq~
    select asmbl_id, sequence, datalength(sequence) from assembly where asmbl_id = $temp[0]~;
# foreach my $orf (@list){

#     my @temp = split('::', $line);

#     $query .= " feat_name like \"$temp[0]\" OR";
# }

# $query =~ s/OR$//;

my $rh = $dbproc->prepare($query)
    || die ("Cannot prepare $query: " . $dbproc->errstr);

$rh->execute()
    || die ("Failed query: " . $dbproc->errstr);

    
	#my $comname;

	while (my $lineref = $rh->fetchrow_arrayref) {

	   my $newline = ">ID" . $$lineref[0] . " " . $$lineref[2] . " bp\n" .  print_sequence2($$lineref[1], 80);

	   $outputfile = "ID" . $$lineref[0];

	       unless (open(FILEOUT, ">$outputfile")) {
		   print STDERR "Can't open file!!\n\n";
		   exit;
	       }

	   print FILEOUT $newline;
	
       }
}

exit;




