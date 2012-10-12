#!/usr/bin/env perl

eval 'exec /usr/bin/perl  -S $0 ${1+"$@"}'
    if 0; # not running under some shell

eval 'exec /usr/bin/perl  -S $0 ${1+"$@"}'
    if 0; # not running under some shell
BEGIN{foreach (@INC) {s/\/usr\/local\/packages/\/local\/platform/}};
use lib (@INC,$ENV{"PERL_MOD_DIR"});
no lib "$ENV{PERL_MOD_DIR}/i686-linux";
no lib ".";


=head1 NAME

prepare_soapdenovo_config.pl - Collates all info and creates config file passed to soapdenovo
                             
=head1 SYNOPSIS

USAGE: ./prepare_soapdenovo_config.pl --max_rd_len=##
				      --avg_ins=##
				      --reverse_seq=0
				      --asm_flags=##
                                      --pair_num_cutoff=##
				      --map_len=##
				      --q1=/path/to/input/prefix.1.fastq
				      --q2=/path/to/input/prefix.2.fastq
				      --q=/path/to/input/prefix.fastq
				      --output_file=/path/to/output_file.txt

=head1 OPTIONS

B<--max_rd_len, -l>
    Maximal read length
    
B<--avg_ins, -i>
    Average insert size
    
B<--reverse_seq, -r>
    Forward/reverse library. 1 if sequence needs to be reversed
    
B<--asm_flags, -a>
    Reads used for contigging and scaffolding.

B<--pair_num_cutoff, -p> 
    Number of mates needed to scaffold across a gap

B<--map_len, -m>
    Minimum length of read mapping to a config

B<--q1>
    Path to input file for read 1 (Fastq format)

B<--q2>
    Path to input file for read 2 (Fastq format)

B<--q>
    Path to input file for singleton reads, if available

B<--output_file, -o>
    Path to output config file

B<--help, -h>
    Print perldocs for this script.
    
=head1 DESCRIPTION

This file exists to prepare the input config that will be fed into the soapdenovo component

=head1 INPUT

 Arguments 

=head1 OUTPUT

 Config file containing all arguments

=head1 CONTACT

    Kemi Abolude
    kabolude@som.umaryland.edu

=cut

use strict;
use warnings;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Pod::Usage;
use File::Basename;
use Ergatis::Logger;
umask(0000);
my $logger;

my %options = &_parse_options();
my $output_file = $options{'output_file'};

open (OUT, "> $output_file") or $logger->logdie("Could not write to $output_file: $!");
print OUT "max_rd_len=$options{'max_rd_len'}\n";
print OUT "[LIB]\n";
print OUT "avg_ins=$options{'avg_ins'}\n";
print OUT "reverse_seq=$options{'reverse_seq'}\n";
print OUT "asm_flags=$options{'asm_flags'}\n";
print OUT "pair_num_cutoff=$options{'pair_num_cutoff'}\n";
print OUT "map_len=$options{'map_len'}\n";
print OUT "q1=$options{'q1'}\n";
print OUT "q2=$options{'q2'}\n";
if ($options{'q'}){
	my $single = `head -n 1 $options{'q'}`;
	print OUT "q=$single\n";
} else{
	print OUT "q=\n";
}
close (OUT);
    
exit(0);

# Parse command-line arguments                                       
sub _parse_options {
    my %opts = ();

    GetOptions(\%opts,
		'max_rd_len|l=i',
		'avg_ins|i=i',
        	'reverse_seq|r:i',
                'asm_flags|a:i',
                'pair_num_cutoff|p:i',
                'map_len|m:i',
                'q1=s',
                'q2=s',
                'q:s',
                'output_file|o=s',
                'help' ) || pod2usage();

    if ($opts{'help'}) {
        pod2usage ( { -exitval => 0, -verbose => 2, -output => \*STDERR } );
    }
    
    my $logfile = Ergatis::Logger::get_default_logfilename();
    my $debug = 4;
    $logger = new Ergatis::Logger( 'LOG_FILE'   =>  $logfile,
                                   'LOG_LEVEL'  =>  $debug );
    $logger = Ergatis::Logger::get_logger();
    

    defined ($opts{'output_file'}) || $logger->logdie("Please specify an output file");
    
    #check that input files exist
    my @input_files = qw( q1 q2 );
    for my $file ( @input_files ) {
        unless (-e $opts{$file}) {
            $logger->logdie("File $opts{$file} does not exist: $!");
        }
        unless (-r $opts{$file}) {
            $logger->logdie("File $opts{$file} is not readable: $!");
        }
        unless (-s $opts{$file}) {
            $logger->logdie("File $opts{$file} is zero-size: $!");
        }

    }
    return %opts
}

