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

USAGE: ./prepare_soapdenovo_config.pl --input_file=/path/to/input/file.list
				      --output_file=/path/to/output_file.txt

=head1 OPTIONS

B<--input_file>
    Path to input file for singleton reads, if available

B<--output_file, -o>
    Path to output config file

B<--help, -h>
    Print perldocs for this script.
    
=head1 DESCRIPTION

This file exists to convert input list files to tab delimited format. Mostly useful for converting clovr tag input to format input acceptable by some components e.g. Bowtie

=head1 INPUT

 List file containing one file path (for single reads) or two file paths (for paired-end reads).
e.g.
/path/to/mate1.fastq
/path/to/mate2.fastq

=head1 OUTPUT

File in tab-delimited format- paired paths are separated by a tab
e.g. 
/path/to/mate1.fastq	/path/to/mate2.fastq


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
my $input_file = $options{'input_file'};
my $output_file = $options{'output_file'};

open (IN, "< $input_file") or $logger->logdie("Could not write to $input_file: $!");
open (OUT, "> $output_file") or $logger->logdie("Could not write to $output_file: $!");
my $tracker = 0;
while (<IN>) {
     chomp;
     if ($tracker == 0) {
	print OUT "$_\t";
     } elsif ($tracker == 1){
	print OUT $_;
     }	
     $tracker++;
}
close (OUT);
close (IN);
    
exit(0);

# Parse command-line arguments                                       
sub _parse_options {
    my %opts = ();

    GetOptions(\%opts,
                'input_file|i=s',
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
    unless (-e $opts{'input_file'}) {
            $logger->logdie("File $opts{'input_file'} does not exist: $!");
    }
    return %opts
}

