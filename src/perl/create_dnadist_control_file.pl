#!/usr/bin/perl

=head1  NAME 

create_dnadist_control_file.pl - Creates a complete dnadist control file from a partial control file.

=head1 SYNOPSIS

create_dnadist_control_file.pl
         --dnadist_seq_file=/path/to/dnadist_aligned_seqs.txt
         --output_control_file=/path/to/complete_control_file.txt
        [--control_file=/path/to/dnadist_partial_control_file.txt
         --log=/path/to/some.log
         --debug=4 ]

=head1 OPTIONS

B<--dnadist_seq_file,-s> 
    path to the file of aligned sequences on which dnadist is to be run.

B<--control_file,-c>
    optional.  path to a dnadist control file that does not contain either: 
    1. an input file name as its first line or 2. a 'Y' as its final line.

B<--output_control_file,-o>
    path to the (complete) output control file.

B<--log,-l> 
    optional.  path to a log file the script should create.  will be overwritten if
    already exists.
    
B<--debug,-d> 
    optional.  the debug level for the logger (an integer)

B<--help,-h> 
    This help message/documentation.

=head1 DESCRIPTION

Creates a complete dnadist control file from a partial control file.  Takes as
its (optional) input a partial control file and then prepends the input file name 
and appends a trailing 'Y' line.

=head1 INPUT

A dnadist sequence file and a partial dnadist control file (i.e., one that does not
set the input sequence file name or end with a trailing 'Y' line.)  See the PHYLIP
dnadist documentation for more information on these file formats.

=head1 OUTPUT

Writes a complete dnadist control file to the file specified by --output_control_file

=head1 CONTACT

    Jonathan Crabtree
    jonathancrabtree@gmail.com

=cut

use strict;
use FileHandle;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Pod::Usage;
use Ergatis::Logger;

## input
my $options = {};
&GetOptions($options,
            "dnadist_seq_file|s=s",
            "control_file|c=s",
            "output_control_file|o=s",
            "log|l=s",
            "debug|d=i",
            "help|h",
            ) || pod2usage();

## display documentation
if( $options->{'help'} ){
    pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDERR} );
}

&check_parameters($options);

## logging
my $logfile = $options->{'log'} || Ergatis::Logger::get_default_logfilename();
my $logger = new Ergatis::Logger('LOG_FILE' => $logfile, 'LOG_LEVEL' => $options->{'debug'});
$logger = $logger->get_logger();

## main program
my $ofh = FileHandle->new();
my($ocf, $sf, $cf) = map {$options->{$_}} qw(output_control_file dnadist_seq_file control_file);
$ofh->open(">$ocf") || die "unable to write to $ocf";
$logger->warn("--dnadist_seq_file=$sf not found") if (!-e $sf);
$ofh->print($sf . "\n");

# echo --control file, if specified
if (defined($cf)) {
    my $ifh = FileHandle->new();
    $ifh->open($cf, 'r') || die "unable to read from $cf";
    while (my $line = <$ifh>) {
        $ofh->print($line);
    }
    $ifh->close();
}

$ofh->print("Y\n");
$ofh->close();
exit(0);

## subroutines
sub check_parameters {
    my $options = shift;

    
    ## make sure required parameters were passed
    my @required = qw(dnadist_seq_file output_control_file);
    for my $option ( @required ) {
        unless  ( defined $options->{$option} ) {
            die "--$option is a required option";
        }
    }

    ## additional parameter checking
    
    ## defaults
}
