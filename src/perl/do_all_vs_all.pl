#!/usr/local/bin/perl

use lib("shared");
use strict;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use English;
use File::Basename;


my %options = ();
my $results = GetOptions (\%options, 'db_file|d=s', 'query_dir|q=s', 'seq_dir|s=s', 'pep_dir|p=s', 'output_dir|o=s',
                                     'asmbl_ids|a=s', 'DEBUG', 'help|h' );

###-------------PROCESSING COMMAND LINE OPTIONS-------------###

my $ASMBL_IDS = $options{'asmbl_ids'};
my $db_file   = $options{'db_file'};
my $query_dir = $options{'query_dir'};
$query_dir =~ s/\/+$//; 
my $seq_dir   = $options{'seq_dir'};
$seq_dir =~ s/\/+$//; 
my $pep_dir   = $options{'pep_dir'};
$pep_dir =~ s/\/+$//; 
my $output_dir     = $options{'output_dir'};
$output_dir =~ s/\/+$//;       

if(!$db_file or !$query_dir or !$seq_dir or !$pep_dir or !$output_dir or !$ASMBL_IDS or exists($options{'help'})) {
    &print_usage();
}

###-------------------------------------------------------###

my @asm_ids;
if(defined($ASMBL_IDS)) {
    @asm_ids = split(/,/, $ASMBL_IDS);
}
foreach my $asmbl_id (@asm_ids) {
    my $error_code = error_check($asmbl_id);
    if($error_code == 0) {
	print "Looks good for asmbl_id $asmbl_id.\nPerforming AllvsAll\n";
	execute_allvsall($asmbl_id);
    }else {
	print "Skipping asmbl_id $asmbl_id, due to error code of $error_code\n";
	exit $error_code;
    }

}


sub execute_allvsall {

    my $asmbl_id = shift;
    

    my $query_fasta = "${asmbl_id}.pep";

    my $fastafile = "$query_dir/$query_fasta";
    my $seqdir    = "$seq_dir/${asmbl_id}";
    my $pepdir    = "$pep_dir";
    my $outputdir = "$output_dir/${asmbl_id}";
	 
	      
    #my $command = "pallvsall -nowait -parameters database=$db_file,fastafile=$fastafile,seqdir=$seqdir,pepdir=$pepdir,outputdir=$outputdir";     
    my $command = "pallvsall -parameters database=$db_file,fastafile=$fastafile,seqdir=$seqdir,pepdir=$pepdir,outputdir=$outputdir";
    print "$command\n";
    qx($command);

    return 1;

}

sub error_check  {

    my $asmbl_id = shift;

    #$org = uc($org);
    my $query_fasta = $asmbl_id.".pep";

    #check for the presence of database file
    if( ! -s $db_file) {
	print "Unable to locate the database file ".basename($db_file)." \n";
	return 2;
    } else {
	qx(setdb $db_file) if(! -e "$db_file.ahd" || ! -e "$db_file.atb" || ! -e "$db_file.bsq");
    }
    
    #check for the presence of query file
    if( ! -d $query_dir) {
	print "Unable to locate the directory holding the query fasta file \"$query_dir\"\n";
	return 3;
    }elsif(! -e "$query_dir/$query_fasta") {
	print "Unable to locate the query file $query_fasta!  Aborting...\n";;
	return 4;
    }
	    
    #check for the presence of *.fsa files
    if( ! -d "$pep_dir/$asmbl_id") {
	print "Unable to find the \"$pep_dir/asmbl_id_${asmbl_id}\" directory. Aborting...\n";
	return 5;
    }elsif(! (my @files = <$pep_dir/$asmbl_id/*.fsa>)) {
	print STDERR "No fsa files found in \"$pep_dir/asmbl_id_$asmbl_id\" directory!  Aborting...\n";
	return 6;
    }    
	    
    #check for the presence of *.seq files
    if( ! -d "$seq_dir/$asmbl_id") {
	print "Unable to find the \"$seq_dir/asmbl_id_${asmbl_id}\" directory. Aborting...\n";
	return 7;
    }elsif(! (my @files = <$seq_dir/$asmbl_id/*.seq>)) {
	print STDERR "No seq files found in \"$seq_dir/asmbl_id_${asmbl_id}\" directory!  Aborting...\n";
	return 8;
    }    
    
    #check for the presence of the output_dir
    if(! -d $output_dir) {
	mkdir $output_dir;
	chmod 0777, $output_dir;
    }
    if(! -d "$output_dir/${asmbl_id}") {
	mkdir "$output_dir/${asmbl_id}";
	chmod 0777, "$output_dir/${asmbl_id}";
    }
        
	


    return 0;

}

sub print_usage {


    print STDERR "SAMPLE USAGE:  do_all_vs_all.pl --db_file something.pep --query_dir q_dir --seq_dir s_dir --pep_dir p_dir --output_dir out_dir -a 1,2\n";
    print STDERR "  --db_file    = database file\n";
    print STDERR "  --query_dir  = dir containing query fasta\n";
    print STDERR "  --seq_dir    = sequence dir\n";
    print STDERR "  --pep_dir    = peptide dir\n";
    print STDERR "  --output_dir = dir to save allvsall results\n";
    print STDERR "  --asmbl_ids (only get sequence belong to particular asmbl_ids)\n";
    print STDERR "  --help = This help message.\n";
    exit 1;

}





