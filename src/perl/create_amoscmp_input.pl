#!/usr/bin/perl

eval 'exec /usr/bin/perl  -S $0 ${1+"$@"}'
    if 0; # not running under some shell
BEGIN{foreach (@INC) {s/\/usr\/local\/packages/\/local\/platform/}};
use lib (@INC,$ENV{"PERL_MOD_DIR"});
no lib "$ENV{PERL_MOD_DIR}/i686-linux";
no lib ".";

=head1 NAME

create_amoscmp_input.pl - takes input from identify_human_contaminates to amoscmp

=head1 SYNOPSIS

USAGE: ./create_amoscmp_input.pl --input=/path/to/for_assembly_output --sff_file=/path/to/sff/file --output_directory=/path/to/output/directory

=head1 OPTIONS

B<--input, -i>
	text file produced by identify_human_contaminates.pl  called something like name.foramoscmp
	
B<--sff_file, -s>
	Raw reads in SFF format. This file should contain the sequences in the FASTA file used as input for nucmer.
	
B<--output_directory, -o>
	Output directory where FASTA and qual files will be written to.
	
B<--log, -l>
	Optional. Log file.
	
B<--debug, -d>
	Optional. Debug level.
	
=head1 DESCRIPTION

=head1 INPUT

Standard nucmer delta output.

=head1 OUTPUT

Both a FASTA file and quality file for the top nucmer alignment hits.

=head1 CONTACT

	Cesar Arze
	carze@som.umaryland.edu
	
=cut							

use strict;
use warnings;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use File::Basename;
use Ergatis::Logger;
use MG::Common;

my %global = MG::Common::initialize();
# TODO: Find a better way to handle hard-coded paths here.

my $logger;
my %options = &parse_options();
my $prefix = (split(/\.\//, $options{'input'}))[-2];

open INPUT, "<$options{input}" or die $!;
open OUTPUT, ">$options{'output_directory'}/create_amoscmp_input.tbl" or die $!;

while (my $line = <INPUT>) {
    chomp($line);
    my ($taxid, $dbfile, $acclist) = split(/\t/, $line);
    open LIST, ">$options{'output_directory'}\/$taxid.list" or die $!;
    print LIST join("\n", split(/,/, $acclist)),"\n";
    close LIST;
    run_system_cmd($global{sffbin}."sfffile -o ".$options{'output_directory'}."/".$taxid.".sff -i  ".$options{'output_directory'}."\/".$taxid.".list ".join(" ",split(/,/,$options{'sff_file'})));
    run_system_cmd("$global{sffbin}sffinfo -s $options{'output_directory'}/$taxid\.sff > $options{'output_directory'}/$taxid\.fsa");
    run_system_cmd("$global{sffbin}sffinfo -q $options{'output_directory'}/$taxid\.sff > $options{'output_directory'}/$taxid\.qual");
    print OUTPUT $options{'output_directory'}."/".$taxid.".fsa","\t",$dbfile,"\n";
}

#########################################################################

## Create the read lists for filtering purposes by sfffile

sub run_system_cmd {
    my $cmd = shift;
    my $res = system($cmd);
        $res = $res >> 8;
        
    unless ($res == 0) {
        $logger->logdie("Could not run $cmd");
    }   
        
    return $res;
}

sub parse_options {
	my %opts = ();
	
	GetOptions(\%opts,
				'input|i=s',
				'sff_file|s=s',
				'output_directory|o=s',
				'log|l=s',
				'debug|d=s') || pod2usage();
				
	if ( $opts{'help'} ) {
	    pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDERR} );
	}
	
    ## Logger needs to be set beforehand.
    my $logfile = $opts{'log'} || Ergatis::Logger::get_default_logfilename();
    $logger = new Ergatis::Logger( 'LOG_FILE' 	=>	$opts{'log'},
       							  'LOG_LEVEL'	=>	$opts{'debug'} );
    $logger = Ergatis::Logger::get_logger();
    
    ## Check some of the parameters...
    defined ($opts{'input'}) || $logger->logdie("Please specify an input file.");
    defined ($opts{'sff_file'}) || $logger->logdie("Please provide an SFF file.");
    defined ($opts{'output_directory'}) || $logger->logdie("Please specify an output directory.");
    defined ($opts{'cutoff'}) || ( $opts{'cutoff'} = "0.75" );
    
    return %opts;			
}
