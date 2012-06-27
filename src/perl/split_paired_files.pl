#!/usr/bin/perl

eval 'exec /usr/bin/perl  -S $0 ${1+"$@"}'
    if 0; # not running under some shell

=head1 NAME

split_paired_files.pl - Splits two ordered files with paths to paired reads into files containing paths to one pair each.

=head1 SYNOPSIS

 USAGE: split_paired_files.pl
    --input_file1=/path/to/paired1.list
    --input_file2=/path/to/paried2.list
    --output_directory=/path/to/output_dir/

=head1 OPTIONS

B<--input_file1,-f>
    Single list file with paths to reads files - pair 1

B<--input_file2,-f>
    Single list file with paths to reads files - pair 2. List must be in same order as input_file1

B<--output_dir,-o>
    The directory to put the split files.

B<--log,-l> 
    Log file

B<--debug> 
    Debug level

B<--help,-h>
    Print this message

=head1  DESCRIPTION

    This script is used to convert two list files containing (ordered) paths to paired read files. 
    The files are split to output new files, each containing the path to a complete pair of read files.

=head1 INPUT
    Two input files are required. Each file contains full paths to the paired reads files.
    
=head1 OUTPUT
    The output is a list of files with each file containing paths to each paired file. The paths are
    separated by a tab as shown below:
    /path/to/pair1.fq     /path/to/pair2.fq

=head1  AUTHOR

    Kemi Abolude
    kabolude@som.umaryland.edu

=cut
    
use strict;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Ergatis::Logger;
use Pod::Usage;

my %options = ();
GetOptions(\%options, 
           'input_file1|i=s',
           'input_file2|j=s',
           'output_dir|o=s',
           'log|l=s',
           'debug=s',
	   'help|h') || pod2usage();
my $input_file1;
my $input_file2;
my $output_dir;

## display documentation
if( $options{'help'} ){
    pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDERR} );
}

## make sure everything passed was peachy
&check_parameters(\%options);

## Setup the logger.  See perldoc for more info on usage
my $logfile = $options{'log'} || Ergatis::Logger::get_default_logfilename();
my $logger = new Ergatis::Logger('LOG_FILE' =>$logfile,
    'LOG_LEVEL'=>$options{'debug'});
$logger = $logger->get_logger();


open(my $f1, "$input_file1");
open(my $f2, "$input_file2");

my @pair1 = <$f1>;
my @pair2 = <$f2>;

my $n=0;

while ($pair1[$n]) {
    chomp ($pair1[$n],$pair2[$n]);
    my $m = $n+1;

    my $filepath = "$output_dir/$m.tab";

    ## take any // out of the filepath
    $filepath =~ s|/+|/|g;

    ## open a new output file
    open (my $ofh, ">$filepath") || $logger->logdie("can't create $filepath\n$!");
    print $ofh "$pair1[$n]\t$pair2[$n]\n";
    $n++;
    close($ofh);
}

close($f1);
close($f2);




sub check_parameters {
    my $options = shift;
    
    ## make sure required arguments were passed
    my @required = qw( input_file1 input_file2 output_dir);
    for my $option ( @required ) {
        unless  ( defined $$options{$option} ) {
            die "--$option is a required option";
        }
    }

    ## make sure the input files exists
    if (! -e $options{'input_file1'}) {
       $logger->logdie("input_file1 $options{'input_file1'} does not exist")
    }
    $input_file1 = $options{'input_file1'};

    if (! -e $options{'input_file2'}) {
	 $logger->logdie("input_file2 $options{'input_file2'} does not exist")
    }
    $input_file2 = $options{'input_file2'};

    ## make sure the output_dir exists
    if (! -e $options{'output_dir'}) {
        $logger->logdie("the output directory passed could not be read or does not exist")
    }
    $output_dir = $options{output_dir};


}



