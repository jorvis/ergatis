#!/usr/local/bin/perl


use strict;
#use Log::Log4perl qw(get_logger);
use DBI;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use English;
use CGC;
#use BsmlBuilder;


my %options = ();
my $results = GetOptions (\%options, 'database|d=s', 'username|u=s', 'password|p=s', 'asmbl_ids|a=s',
                                     'simple_header|s', 'output_dir|f=s', 'DEBUG', 'help|h', 'each_file|e' );


###-------------PROCESSING COMMAND LINE OPTIONS-------------###

my $database  = $options{'database'};
my $user      = $options{'username'};
my $password  = $options{'password'};
my $ASMBL_IDS = $options{'asmbl_ids'};
my $simple_header = $options{'simple_header'};
my $is_public = $options{'ispublic'} || 1;
my $each_file = $options{'each_file'} || 0;
my $output_dir = $options{'output_dir'};
my $QUERYPRINT;
my $DEBUG = $options{'DEBUG'} || 0;

#Log::Log4perl->init("log.conf");
#my $logger = get_logger();



if(!$database or !$user or !$password or exists($options{'help'})) {
    #$logger->fatal("Not all of the required options have been defined.  Exiting...");
    &print_usage();
}

###-------------------------------------------------------###

#Make new CGC object#
my $CGC = new CGC( user       => $user,
                   password   => $password,
                   db         => $database,
                   debug      => $DEBUG,
		   queryprint => $QUERYPRINT
		 );

=hello
my $output_dir = $ENV{'PNEUMO_DB_DIR'} if(!$output_dir);
$output_dir =~ s/\/+$//;       #remove terminating '/'s
if(! -d $output_dir ) {
    mkdir $output_dir;
}
chmod 0777, $output_dir;
=cut

my @asm_ids;
if(defined($ASMBL_IDS)) {
    @asm_ids = split(/,/, $ASMBL_IDS);
}else {
    @asm_ids = get_asmbl_id();
}


foreach my $asmbl_id (@asm_ids) {
    #my $final_output_dir;
    my $asmbl_sequence = get_asmbl_seq($asmbl_id);
    if(!defined($asmbl_sequence)) {
	print STDERR "No asmbl_sequence can be found for asmbl_id $asmbl_id\n";
	next;
    } else {
	my $fasta_header = "PNEUMO_asmbl_$asmbl_id";
	my $fastaout = &fasta_out($fasta_header, $asmbl_sequence);
	print $fastaout;
    }
}


sub fasta_out {
#This subroutine takes a sequence name and its sequence and
#outputs a correctly formatted single fasta entry (including newlines).

    my $seq_name = shift;
    my $seq = shift;

    my $fasta=">"."$seq_name"."\n";
    for(my $i=0; $i < length($seq); $i+=60){
	my $seq_fragment = substr($seq, $i, 60);
	$fasta .= "$seq_fragment"."\n";
    }
    return $fasta;

}


sub get_asmbl_seq {
#This subroutine returns the sequence of the assembly given asmbl_id
#The returned data is a string representing the sequence of an asmbl_id

    my $asmbl_id = shift;

    my $result = $CGC->fetch_assembly_sequence_of($asmbl_id);
    my $asmbl_seq = $result->{0}->{'sequence'};


    return $asmbl_seq;

}


sub get_asmbl_id {


    my $result = $CGC->fetch_all_asmbl_ids();
    
    my @ASMBL_IDS;
    for(my $i=0; $i<$result->{'count'}; $i++) {
	my $asmbl_id = $result->{$i}->{'asmbl_id'};
	push(@ASMBL_IDS, $asmbl_id);
    }

    return @ASMBL_IDS;

}
