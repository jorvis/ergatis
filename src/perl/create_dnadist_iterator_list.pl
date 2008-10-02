#!/usr/bin/perl

=head1  NAME 

create_dnadist_iterator_list.pl - Create ergatis/workflow iterator for distributed dnadist job.

=head1 SYNOPSIS

create_dnadist_iterator_list.pl 
         --dnadist_seq_file=/path/to/dnadist_aligned_seqs.txt
         --output_iterator_list=/path/to/output_iterator_list.txt
        [--group_size=1000
         --log=/path/to/some.log
         --debug=4 ]

=head1 OPTIONS

B<--dnadist_seq_file,-s> 
    path to the file of aligned sequences on which dnadist is to be run. (used only
    to determine how many sequences there are.)

B<--output_iterator_list,-o> 
    path to the ergatis/workflow iterator list file to create.

B<--group_size,-g> 
    optional.  specifies the maximum number of (pairwise) sequence distance calculations
    to perform in each iterator group.  Defaults to 1000.  Note that for an input of N
    sequences, (N(N-1)) / 2 pairwise comparisons will be performed.

B<--log,-l> 
    optional.  path to a log file the script should create.  will be overwritten if
    already exists.
    
B<--debug,-d> 
    optional.  the debug level for the logger (an integer)

B<--help,-h> 
    This help message/documentation.

=head1 DESCRIPTION

Creates an ergatis/workflow iterator list file for a distributed dnadist job.

=head1 INPUT

A sequence file in the format accepted by dnadist. The script currently looks only at
the first line of the file, which is expected to have two whitespace-delimited fields:

100 7682

The first field is the number of sequences in the file and the second field is the 
number of characters (i.e., the sequence length) in each aligned sequence.  The number
of sequences is used to determine the size of the output matrix and divide up the 
cells of the matrix amongst the dnadist_worker processes.

=head1 OUTPUT

The output ergatis iterator list will be written to the file specified by --output_iterator_list

=head1 CONTACT

    Jonathan Crabtree
    jonathancrabtree@gmail.com

=cut

use strict;
use FileHandle;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Pod::Usage;
use Ergatis::Logger;

## globals
my $DEFAULT_GROUP_SIZE = 1000;

## input
my $options = {};
&GetOptions($options,
            "dnadist_seq_file|s=s",
            "output_iterator_list|o=s",
            "group_size|g=i",
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

# parse first line of input file to determine number of input sequences
my $nseqs = undef;
my $nchars = undef;

my($sf, $il, $gs) = map {$options->{$_}} qw(dnadist_seq_file output_iterator_list group_size);

my $ofh = FileHandle->new();
$ofh->open(">$il") || die "unable to write to $il";

my $ifh = FileHandle->new();
$ifh->open($sf, 'r') || die "unable to read from $sf";
my $line = <$ifh>;
if ($line =~ /^\s*(\d+)\s+(\d+)\s*$/) {
    $nseqs = $1;
    $nchars = $2;
}
$ifh->close();

die "unable to determine number of sequences in input file ($sf)" if (!defined($nseqs));

# calculate number of comparisons
my $ncomparisons = ($nseqs * ($nseqs-1)) / 2;

# ergatis iterator file header
$ofh->print(join("\t", '$;I_START_CELL$;', '$;I_END_CELL$;') . "\n");

# divide $ncomparisons into groups of size <= $gs
for (my $first = 0;$first < $ncomparisons;$first += $gs) {
    my $last = $first + $gs - 1;
    $last = $ncomparisons - 1 if ($last > ($ncomparisons-1));
    $ofh->print(join("\t", $first, $last) . "\n");
}

$ofh->close();
exit(0);

## subroutines
sub check_parameters {
    my $options = shift;
    
    ## make sure required parameters were passed
    my @required = qw(dnadist_seq_file output_iterator_list);
    for my $option ( @required ) {
        unless  ( defined $options->{$option} ) {
            die "--$option is a required option";
        }
    }

    ## additional parameter checking
    my $sf =  $options->{'dnadist_seq_file'};
    if (!-e $sf) {
        die "--dnadist_seq_file=$sf does not exist";
    }

    my $il =  $options->{'output_iterator_list'};
    if (-e $il) {
        die "--output_iterator_list=$il already exists";
    }
    
    ## defaults
    if (!defined($options->{'group_size'} || ($options->{'group_size'} <= 0))) {
        $options->{'group_size'} = $DEFAULT_GROUP_SIZE;
    }
}
