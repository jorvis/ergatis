#!/usr/bion/perl

=head1 NAME

./replace_amoscmp_fasta_seqids.pl - takes FASTA output from amoscmp and replaces the sequence identifiers with unique counterparts

=head1 SYNOPSIS

./replace_amoscmp_fasta_seqids.pl
        --fasta_input_file=/path/to/fasta/file
        --log=/path/to/log/file
        --debug=/debug/level
        --help

=head1 PARAMETERS

B<--fasta_input_file, -i>
    An input FASTA file produced by amoscmp. The basename of this input file will be utilized
    when generating unique sequence identifiers
    
B<--log, -l>
    Optional. Log file

B<--debug, -d>
    Optional. Debug level
    
B<--help>
    Prints out script documentation.
    
=head1 DESCRIPTION

This script aims to replace the default sequence ids produce by amoscmp in its FASTA output from
generic numbers to one containing the basename of the input file. This should help prevent any
duplication errors that may be encountered in down stream analysis.

    Input File: ecol_contig3.fsa
    > 1     -->     > ecoli_contig3_1
    > 2     -->     > ecoli_contig3_2
    > 3     -->     > ecoli_contig3_3
  
=head1 INPUT

The FASTA output file created from an amoscmp assembly. 

=head1 OUTPUT
   
A FASTA file containing sequence ID's less generic and less prone to cause duplication errors
in downstream analysis.

=head1 CONTACT

    Cesar Arze
    carze@som.umaryland.edu
    
=cut

use strict;
use File::Basename;
use Ergatis::Logger;
use Tie::File;
use Pod::Usage;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);

## GLOBALS
my $logger;

my %options = &parse_options();
my $fasta_file = $options{'fasta_input_file'};

## The basename of the FASTA file will be used as one part 
## of the new sequence identifier
my $filename = fileparse($fasta_file, '\.[^\.]*'); 

## Tie FASTA file to a hash so we can edit in place
my @fasta = ();
tie @fasta, 'Tie::File', $fasta_file or $logger->logdie("Could not open FASTA file $fasta_file: $!");
foreach my $line (@fasta) {
    ## If we are dealing with a header we want to replace it with
    ## our new sequence identifier
    if ($line =~ /^>(\d+)/) {
        my $new_seq_id = $filename . "_" . $1;
        $line = ">$new_seq_id\n";
    }
}

untie @fasta;

###############################################################################
#####                          SUBROUTINES                                #####
###############################################################################

sub parse_options {
    my %opts = ();
    GetOptions(\%opts,
                'fasta_input_file|f=s',
                'log|l:s',
                'debug|d:s',
                'help=s') || pod2usage();
                
	&pod2usage( {-exitval => 1, -verbose => 1, -output => \*STDOUT}) if ($opts{'help'} );                
    
    ## Initialize logger
    my $log_file = $options{'log'} || Ergatis::Logger->get_default_logfilename();
    my $debug_lvl = $options{'debug'} ||= 4;
    $logger = new Ergatis::Logger( 'LOG_FILE' => $log_file, 'LOG_LEVEL' => $debug_lvl );
    $logger= Ergatis::Logger->get_logger();
    
    defined ($opts{'fasta_input_file'}) || $logger->logdie("Please specify a valid input FASTA file");
    
    return %opts;
}