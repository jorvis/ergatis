#!/usr/local/bin/perl

use strict;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use File::Basename;


my %options = ();
my $results = GetOptions (\%options, 'db_file|d=s', 'query_fasta|q=s', 'output_dir|o=s', 'help|h' );

###-------------PROCESSING COMMAND LINE OPTIONS-------------###

#my $ASMBL_IDS = $options{'asmbl_ids'};
my $db_file   = $options{'db_file'};
my $query_fasta = $options{'query_fasta'};
my $output_dir     = $options{'output_dir'};
$output_dir =~ s/\/+$//;       #remove terminating '/'s

if(!$db_file or !$query_fasta or !$output_dir or exists($options{'help'})) {
    print STDERR "Require db_file, query_fasta, outputdir\n";
    exit 10;
}

###-------------------------------------------------------###


#check for the presence of database file and query file
if( ! -s $db_file) {
    print "Unable to locate the database file $db_file. Aborting...\n";
    exit 2;
} else {
    qx(setdb $db_file) if(! -e "$db_file.ahd" || ! -e "$db_file.atb" || ! -e "$db_file.bsq");
}

if( ! -s $query_fasta) {
    print "Unable to locate the query fasta file $query_fasta. Aborting... \n";
    exit 3;
}
    
if(! -d $output_dir) {
    mkdir $output_dir;
} else {
    unlink glob("$output_dir/*");  #delete old output files if any....
}
chmod 0777, $output_dir;



my $command = "pblastp -parameters database=$db_file,fastafile=$query_fasta,report=both,outputdir=$output_dir";
qx($command);




