#!/usr/local/bin/perl

use lib (@INC,$ENV{"PERL_MOD_DIR"});
no lib "$ENV{PERL_MOD_DIR}/i686-linux";
no lib ".";

=head1  NAME 

die_if_empty.pl  - takes a list of input files and exits with an error if any are empty

=head1 SYNOPSIS

USAGE:  die_if_empty.pl -i /path/to/input.list [--not]

=head1 OPTIONS

B<--input_list,-i>
    The list(s) of files to check for emptiness [can be two or more comma delimited]

B<--input_file,-f>
    The file to check
    
B<--not>
    Optionally exit with an error if the files are NOT empty.

B<--log,-l> 
    Log file

B<--debug,-d>
    Debug level
    
B<--help,-h>
    This help message

=head1   DESCRIPTION

Allows checks for emptiness or non-emptiness of a set of files.

=head1  CONTACT
    Brett Whitty
    bwhitty@tigr.org

=cut

use strict;
use warnings;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use Pod::Usage;
BEGIN {
use Ergatis;:Logger;
}

my %options = ();
my $results = GetOptions (\%options, 
                          'input_list|i:s',
                          'input_file|f:s',
                          'log|l=s',
                          'debug|d=s',
                          'not',
                          'help|h') || pod2usage();

my $logfile = $options{'log'} || Ergatis;:Logger::get_default_logfilename();
my $logger = new Ergatis;:Logger('LOG_FILE'=>$logfile,
                                  'LOG_LEVEL'=>$options{'debug'});
$logger = $logger->get_logger();

my @files;
my $empty_count=0;
my $nonempty_count=0;

## display documentation
if( $options{'help'} || scalar keys(%options) == 0 ){
    pod2usage( {-exitval=>0, -verbose => 2, -output => \*STDOUT} );
}

if ($options{'input_list'}) {
    my @lists = split(",", $options{'input_list'});
    foreach my $list(@lists) {
        unless(-e $list) {
            $logger->logdie("provided input list '$list' doesn't exist");
        }

        open(LIST, $list) || $logger->logdie("couldn't open input list '$list' for reading");

        while (<LIST>) {
            chomp;
            if ($logger->is_debug) { $logger->debug("adding file '$_' for checking");}
            push (@files, $_);
        }
    }
}

if ($options{'input_file'}) {
    unless(-e $options{'input_file'}) {
        $logger->logdie("provided input list $options{input_file} doesn't exist");
    }
    if ($logger->is_debug) { $logger->debug("adding file '$options{input_file}' for checking");}
    push (@files, $options{'input_file'});
}

foreach my $file(@files) {
    if (-z $file) {
        $empty_count++;
        if ($logger->is_debug) { $logger->debug("file '$file' is empty");}
    } else {
        $nonempty_count++;
        if ($logger->is_debug) { $logger->debug("file '$file' is not empty");}
    }
}

if      ( $options{'not'}   &&    $nonempty_count   > 0) {
    $logger->logdie("one or more of the specified files were not empty");
} elsif (!$options{'not'}   &&    $empty_count      > 0) {
    $logger->logdie("one or more of the specified files were empty");
}
