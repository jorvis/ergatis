#!/usr/local/bin/perl

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
    execute_allvsall($asmbl_id);
}

exit(0);

sub execute_allvsall {

    my $asmbl_id = shift;
    

    my $query_fasta = "${asmbl_id}.pep";

    my $fastafile = "$query_dir/$query_fasta";
    my $seqdir    = "$seq_dir/${asmbl_id}";
    my $pepdir    = "$pep_dir";
    my $outputdir = "$output_dir/${asmbl_id}";
	 
    if(! -d $outputdir) {
	mkdir $outputdir;
	chmod 0777, $outputdir;
    }
	      
    if( !( -s $fastafile ))
    {
	# The assembly does not have any genes
	exit(0);
    }
     
    my $command = "pallvsall -parameters database=$db_file,fastafile=$fastafile,seqdir=$seqdir,pepdir=$pepdir,outputdir=$outputdir";
    my $status = system( $command );

    my $exit_value = $status >> 8;
    my $signal_num = $status & 127;
    my $dumped_core = $status & 128;

    if( !($exit_value == 0) )
    {
	exit( $exit_value );
    }

    return $exit_value;
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





