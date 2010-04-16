#!/usr/bin/perl

eval 'exec /usr/bin/perl  -S $0 ${1+"$@"}'
    if 0; # not running under some shell
BEGIN{foreach (@INC) {s/\/usr\/local\/packages/\/local\/platform/}};
use lib (@INC,$ENV{"PERL_MOD_DIR"});
no lib "$ENV{PERL_MOD_DIR}/i686-linux";
no lib ".";

=head1 NAME

sfffile.pl - Wrapper script for sfffile to facilitate its execuction via Ergatis

=head1 SYNOPSIS

./sfffile.pl 
        --sfffile_exe=/path/to/sfffile/executable
        --raw_sff_files=/path/to/sff/files
        --command_line_opts=/sfffile/command/line/options
        --accession_list=/path/to/accession/list
        --output_file=/path/to/target/output/file
        --log=/path/to/log/file
        --debug=/debug/level
        --help

=head1 PARAMETERS

B<--sfffile_exe, -e>
    Path to the sfffile executable

B<--raw_sff_files, -s> 
    Either a single SFF file or a comma-delimited list of SFF files.
    
B<--command_line_opts, -c>
    Any command line options that should be passed to the sfffile executable
    
B<--accession_list, -a>
    A list of accessions that should be used to generate the new filtered SFF file
    
B<--output_file, -o>
    Path to the target output SFF file.
    
B<--log, -l>
    Optional. Log file
    
B<--debug, -d>
    Optional. Debug level
    
B<--help>
    Print out documentation for this wrapper script
    
=head1 DESCRIPTION

A wrapper script that handles executing the sfffile program. Specifically it can handle a comma-delimited list of SFF files
provided by the user.

=head1 INPUT

A list of accessions to filter with and a comma-delimited list of SFF files

=head1 OUTPUT

A new SFF containing only the reads corresponding to the accession list

=head1 CONTACT

    Cesar Arze
    carze@som.umaryland.edu
    
=cut

use strict;
use warnings;
use Pod::Usage;
use Ergatis::Logger;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
my $logger;

my %options = ();
my $results = GetOptions(\%options,
                          'raw_sff_files|s=s',
                          'accession_list|a=s',
                          'sfffile_exe|e=s',
                          'output_file|o=s',                          
                          'command_line_opts|c:s',
                          'log|l=s',
                          'debug|d=s',
                          'help=s') || pod2usage();
                          
&pod2usage( {-exitval => 1, -verbose => 1, -output => \*STDOUT}) if ($options{'help'} );                
    
## Initialize logger
my $log_file = $options{'log'} || Ergatis::Logger->get_default_logfilename();
my $debug_lvl = $options{'debug'} ||= 4;
$logger = new Ergatis::Logger( 'LOG_FILE' => $log_file, 'LOG_LEVEL' => $debug_lvl );
$logger= Ergatis::Logger->get_logger();

## Check if the sfffile executable, SFF files, and accession list were provided...
defined ($options{'sfffile_exe'}) || $logger->logdie("Pleasep provide the path to the sfffile executable.");
defined ($options{'raw_sff_files'}) || $logger->logdie("Please provide one or more SFF files for processing.");
defined ($options{'accession_list'}) || $logger->logdie("Please provide a list of accessions to filter with.");
defined ($options{'output_file'}) || $logger->logdie("Please provide a valid output file.");

my $sfffile_exec = $options{'sfffile_exe'};
my $sff_files = $options{'raw_sff_files'};
my $accessions = $options{'accession_list'};
my $output_file = $options{'output_file'};
my $command_line_opts;

if (defined($options{'command_line_opts'})) {
    $command_line_opts = $options{'command_line_opts'}
} else {
    $command_line_opts = "";
}

## If a comma-delimited list of SFF files are passed in we want to change it to a space-delimited list
$sff_files =~ s/\,/ /g;

my $cmd = "$sfffile_exec -o $output_file -i $accessions $command_line_opts $sff_files";

my $res = system($cmd);
$res = $res >> 8;
       
unless ($res == 0) {
    $logger->logdie("Could not run $cmd");
}
                          

    
                                     