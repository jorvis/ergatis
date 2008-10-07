#!/usr/bin/perl

=head1  NAME 

join_multifasta.pl  - joins a set of fasta input files into one multifasta file

=head1 SYNOPSIS

USAGE:  join_multifasta.pl -i /path/to/input.list -f /path/to/input_file.fsa -D /path/to/files -o /path/to/output/file.fsa --compress 0

=head1 OPTIONS

B<--input_file_list,-i>
    The list(s) of files to join in a multifasta file [can be two or more comma delimited]

B<--input_file,-f>
    The file to check

B<--input_dir,-D>
    The input directory

B<--input_directory_extension,-s>
    The input file suffix

B<--output_file,-o>
    The output file name

B<--compress>
    Compress the output files (1 = compress / 0 = no compression)
    
B<--log,-l> 
    Log file

B<--debug,-d>
    Debug level
    
B<--help,-h>
    This help message

=head1   DESCRIPTION

Joins a set of fasta files into a single multifasta file.

=head1  CONTACT
    Brett Whitty
    bwhitty@jcvi.org

=cut

use strict;
use warnings;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use Pod::Usage;
use File::Find;
use Ergatis::Logger;

## required for proper operation on NFS
##  see SF.net bug 2142533 - https://sourceforge.net/tracker2/?func=detail&aid=2142533&group_id=148765&atid=772583
$File::Find::dont_use_nlink = 1;

my %options = ();
my $results = GetOptions (\%options, 
                          'input_file_list|i:s',
                          'input_file|f:s',
                          'input_directory|D:s',
                          'input_directory_extension|s=s',
                          'output_file|o=s',
                          'compress=i',
                          'log|l=s',
                          'debug|d=s',
                          'help|h') || pod2usage();

my $logfile = $options{'log'} || Ergatis::Logger::get_default_logfilename();
my $logger = new Ergatis::Logger('LOG_FILE'=>$logfile,
                                  'LOG_LEVEL'=>$options{'debug'});
$logger = $logger->get_logger();

my @files;

## display documentation
if( $options{'help'} || scalar keys(%options) == 0 ){
    pod2usage( {-exitval=>0, -verbose => 2, -output => \*STDOUT} );
}

my @infiles = get_input_files();

my $out_fh = open_file_write($options{'output_file'}, $options{'compress'});

foreach my $file(@infiles) {
    my $in_fh = open_file_read($file);
    while (<$in_fh>) {
        print $out_fh $_;
    }
    close $in_fh;
}
close $out_fh;

## opens a filehandle for reading
sub open_file_read {
    my ($filename) = @_;

    my $fh;

    if (! -e $filename && -e $filename.'.gz') {
        $filename .= '.gz';
    }
    if ($filename =~ /\.gz$|\.gzip$/) {
        open ($fh, "<:gzip", $filename) || die "Couldn't open '$filename' for reading: $!";
    } else {
        open ($fh, "<$filename") || die "Couldn't open '$filename' for reading: $!";
    }

    return $fh
}

## opens a filehandle for writing
sub open_file_write {
    my ($filename, $gzip_mode) = @_;

    my $fh;

    if ($gzip_mode) {
        open ($fh, ">:gzip", $filename.'.gz') || die "Couldn't open '$filename' for writing: $!";
    } else {
        open ($fh, ">$filename") || die "Couldn't open '$filename' for writing: $!";
    }

    return $fh
}

# get a list of input files
sub get_input_files {
    my @infiles;

    if ($options{'input_file_list'}) {
        my @lists = split(",", $options{'input_file_list'});
        foreach my $list(@lists) {
            unless(-e $list) {
                $logger->logdie("provided input list '$list' doesn't exist: $!");
            }

            open(LIST, $list) || $logger->logdie("couldn't open input list '$list' for reading");

            while (<LIST>) {
                chomp;
                if ($logger->is_debug) { $logger->debug("adding file '$_'");}
                push (@infiles, $_);
            }
        }
    }

    if ($options{'input_directory'}) {
        unless (-d $options{'input_directory'}) {
            die "specified input dir '$options{input_directory}' is not a directory";
        }
        if (!$options{'input_directory_extension'}) {
            $options{'input_directory_extension'} = '.fsa';
        }
        find({wanted => sub{ if(/$options{input_directory_extension}(.gz)?$/) {push(@infiles,$_)}}, no_chdir => 1}, ($options{'input_directory'}));
    }

    if ($options{'input_file'}) {
            if (-e $options{'input_file'}) {
                push (@infiles,$options{'input_file'});
            } else {
                print STDERR "ERROR: specified input file '$options{input_file}' doesn't exist.";
            }
    }

    return @infiles;
}
