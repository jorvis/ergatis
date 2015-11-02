#!/usr/bin/perl

eval 'exec /usr/bin/perl  -S $0 ${1+"$@"}'
    if 0; # not running under some shell
BEGIN{foreach (@INC) {s/\/usr\/local\/packages/\/local\/platform/}};
use lib (@INC,$ENV{"PERL_MOD_DIR"});
no lib "$ENV{PERL_MOD_DIR}/i686-linux";
no lib ".";

=head1 NAME

fastq_to_CA.pl - runs the fastqToCA utility found in the wgs-assembler software package. This script acts as a wrapper to ensure that input and 
                 command-line parameters are handled correctly while running this utility through the Ergatis framework. 

=head1 SYNOPSIS

USAGE: fastq_to_CA.pl
            --fastq_to_CA_exe=/path/to/fastqToCA/binary
            --fastq_file=/path/to/fastq/file
            --fastq_file_list=/path/to/fastq/file/list
            --library=<library name>
            --type=<fastq type>
            --insert_size=<mate insert size>
            --reads_orientation=<paired end reads orientation>
           [--log=/path/to/log/file
            --debug=<debug level>]
            
=head1 OPTIONS

B<--fastq_file, -f> 
    A single fastq input file.

B<--fastq_file_list, -fl>
    A tab delimited list of fastq mated reads.
    
B<--fastq_to_CA_exe, -exe>
    Path to the fastqToCA executable.
    
B<--library, -lib>
    The UID of the library these reads are added to.
    
B<--type, -t>
    fastq type; can be one of the three options:
        * sanger - QV's are PHRED, offset=33 '!', NCBI SRA data.
        * solexa - QV's are PHRED, offset=33 '!', NCBI SRA data.
        * illumina - QV's are PHRED, offset=64 '@', Illumina reads from version 1.3 on.

B<--insert_size, -is>
    Mates are on average i +- d bp apart. Must provide the i and d numeric parameters

B<--reads_orientation, -r>
    Orientation of paired-end reads. If 5'-3' <-> 3'-5' set parameter to 'innie.' 
    If 3'-5' <-> 5'-3' set parameter to 'outtie'

B<--output, -o>
    Output path where fragment files will be placed.
    
B<--log, -l>
    Log file

B<--debug>
    Debug level
            
=head1 DESCRIPTION

Wrapper for the fastqtoCA executable. Ensures that input and parameters are correctly formatted for use in the Ergatis framework.

=head1 INPUT

Either a single or list of fastq files that will comprise a single assembly.

=head1 OUTPUT 

Fragment files that are ready to be supplied to the celera-assembler

=head1 CONTACT
    
    Cesar Arze
    carze@som.umaryland.edu

=cut                             

use warnings;
use strict;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Ergatis::Logger;

my %options = ();
my $results = GetOptions (\%options,
                          'fastq_file|f=s',
                          'fastq_file_list|fl=s',
                          'fastq_to_CA_exe|exe=s',
                          'library|lib=s',
                          'type|t=s',,
                          'insert_size|is=s',
                          'reads_orientation|r=s',
                          'output|o=s',
                          'log|l=s',
                          'debug=s',
                          'help|h') || pod2usage();

my $logfile = $options{'log'} || Ergatis::Logger::get_default_logfilename();
my $logger = new Ergatis::Logger( 'LOG_FILE'    =>  $logfile,
                                  'LOG_LEVEL'   =>  $options{'debug'} );
$logger = Ergatis::Logger::get_logger();
                          
my $FASTQ_TO_CA_EXE = $options{'fastq_to_CA_exe'};                        
                          
my $fastq_file = $options{'fastq_file'};
my $fastq_file_list = $options{'fastq_file_list'};
my $library = $options{'library'};
my $fastq_type = $options{'type'};
my $orientation = $options{'reads_orientation'};
my $insert_size = $options{'insert_size'};
my $output = $options{'output'};

my @exe_params = &setup_exe_params($output, $library, $fastq_type, $orientation, $insert_size);                          
my $params = join(' ', @exe_params);                          
my @input;      
                          
if ($options{'fastq_file_list'}) {
    @input = &parse_input_file_list($fastq_file_list);
} elsif ($options{'fastq_file'}) {
    push (@input, $fastq_file);
}

foreach my $input_file (@input) {
    my $cmd = "$FASTQ_TO_CA_EXE $params -fastq $input_file > $output";
    run_system_cmd($cmd);
}

#####################################################
#                                                   #
#                   SUBROUTINES                     #
#                                                   #
#####################################################

sub run_system_cmd {
    my $cmd = shift;
    my $res = system($cmd);
    $res = $res >> 8;
    
    unless ($res == 0) {
        $logger->logdie("Could not run $cmd");
    }
}

sub setup_exe_params {
    my ($out, $lib, $type, $orientation, $ins_size) = @_;
    my @params;
    
    unless ($out) {
        $logger->logdie("Please provide a valid output");
    }
    
    if ($lib) {
    	push (@params, "-libraryname $lib");
    } else {
    	$logger->logdie("Please provide a valid library name");
    }
    
    if ($type) {
        push (@params, "-type $type"); 
    }
    
    if ($orientation) {
        push (@params, "-$orientation");
    }
    
    if ($ins_size) {
        push (@params, "-insertsize $ins_size");
    }
    
    return @params;
}

sub parse_input_file_list {
    my $file = shift;
    my @fastq_mated_reads;
    
    open (FILELIST, $file) or $logger->logdie("Could not open fastq file list $file: $!");
    while (my $line = <FILELIST>) {
        chomp($line);
        my ($mate1, $mate2) = split(/\t/, $line);     
        verify_file($mate1);
        verify_file($mate2);
        
        push (@fastq_mated_reads, "$mate1,$mate2");
    }
    
    return @fastq_mated_reads;
}                 

sub verify_file {
    my $file = shift;
    
    if ( ( defined($file) ) && (-e $file) && (-r $file) && (-s $file) ) {
        return 1;
    } else {
        
        if ( !defined($file) ) {
            $logger->logdie("FASTQ file list $file was not defined");
        }
        
        if ( !-e $file ) {
            $logger->logdie("FASTQ file list $file does not exist");
        }
        
        if ( !-r $file ) {
            $logger->logdie("FASTQ file list $file does not have read permissions");
        }
        
        if ( !-s $file ) {
            $logger->logdie("FASTQ file list $file has zero content");
        }
    }
}
