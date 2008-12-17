#!/usr/bin/perl

=head1 NAME 

concat_and_sort_dnadist_worker_output.pl - Concatenates and sorts intermediate output from dnadist_worker processes.

=head1 SYNOPSIS

concat_and_sort_dnadist_worker_output.pl
         --worker_output_dir=/path/to/dnadist_worker_output
         --output_file=/path/to/complete_control_file.txt
        [--worker_output_regex='^\S+\.distances$'
         --tmp_dir=/tmp
         --worker_output_is_sorted
         --delete_worker_files
         --delete_tmp_files
         --log=/path/to/some.log
         --debug=4 ]

=head1 OPTIONS

B<--worker_output_dir,-w> 
    path to the (parent) directory in which all of the dnadist_worker output
    was written

B<--output_file,-o>
    file where the concatenated and sorted output should be written.

B<--worker_output_regex,-r>
    optional.  regular expression used to identify dnadist_worker output files.
    defaults to \S+\.distances

B<--tmp_dir,-t>
    optional.  directory in which to place temporary files.  defaults to /tmp

B<--worker_output_is_sorted,-a>
    optional.  whether the worker output files are already sorted.

B<--delete_worker_files,-x>
    optional.  whether to delete worker output files once they are no longer needed.

B<--delete_tmp_files,-y>
    optional.  whether to delete temporary files once they are no longer needed.

B<--log,-l> 
    optional.  path to a log file the script should create.  will be overwritten if
    already exists.
    
B<--debug,-d> 
    optional.  the debug level for the logger (an integer)

B<--help,-h> 
    This help message/documentation.

=head1 DESCRIPTION

Concatenates and sorts intermediate output from dnadist_worker processes.
Modified version of the script that assumes the output will be in lower triangular 
form.

=head1 INPUT

A set of files containing whitespace-delimited output of the following form:

1 0.244838
100 0.244838
2 0.056921
200 0.056921

The first value on each line is the 0-based offset of a cell in the output
matrix and the second value is the distance value in that cell.  The output 
must be redundant (i.e. not restricted to the lower triangle of the matrix) 
if the final output will be the full matrix.  If only a lower triangular 
matrix will be output (by dnadist_combiner) then the intermediate input can
omit the duplicate values, like so:

1 0.244838
2 0.056921

=head1 OUTPUT

Writes a file containing all the lines from each of the input files, 
sorted by ascending matrix cell number (the first integer on each line.)

=head1 CONTACT

    Jonathan Crabtree
    jonathancrabtree@gmail.com

=cut

use strict;
use FileHandle;
use File::Basename;
use File::Find;
use File::Spec;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Ergatis::Logger;
use Pod::Usage;

## globals
my $DEFAULT_TMP_DIR = '/tmp';
my $DEFAULT_WORKER_OUTPUT_REGEX = '^\S+\.distances$';

## input
my $options = {};
&GetOptions($options,
            'worker_output_dir|w=s',
            'worker_output_regex|r=s',
            'output_file|o=s',
            'tmp_dir|t=s',
            'worker_output_is_sorted|a',
            'delete_worker_files|x',
            'delete_tmp_files|y',
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

# use File::Find to locate dnadist_worker output files
$File::Find::dont_use_nlink = 1;

my $wd = $options->{'worker_output_dir'};
my @paths = ();
my $regex = $options->{'worker_output_regex'};

my $callback = sub {
    my $filename = $File::Find::name;
    $logger->debug("File::Find checking $filename");
    if ((-f $filename) && ($filename =~ /$regex/)) {
        $logger->debug("adding $filename to file list");
        push(@paths, $filename);
    }
};
&find($callback, $wd);
my $nf = scalar(@paths);
$logger->info("found $nf output file(s) in $wd");
exit(1) if ($nf == 0);

# sort each file individually
my $sorted_files = [];
if ($options->{'worker_output_is_sorted'}) {
    $logger->info("skipping initial sort step due to --worker_output_is_sorted flag");
}

foreach my $path (@paths) {
    if ($options->{'worker_output_is_sorted'}) {
        push(@$sorted_files, $path);
    } else {
        $logger->info("sorting $path\n");
        my($filename, $directories, $suffix) = fileparse($path);
        my $sorted_path = File::Spec->catfile($options->{'tmp_dir'}, $filename . '.sorted');
        &sort_file($path, $sorted_path);
        push(@$sorted_files, $sorted_path);
    }
}

# concatenate sorted files in the correct order
my $file2start = {};
foreach my $file (@$sorted_files) {
    if ($file =~ /(\d+)\-\d+\.distances$/) {
        my $start = $1;
        $file2start->{$file} = $start;
    } else {
        die "error parsing start cell from $file";
    }
}

my @files = sort { $file2start->{$a} <=> $file2start->{$b} } @$sorted_files;
my $first = 1;
my $of = $options->{'output_file'};

foreach my $file (@files) {
    my $cmd = undef;

    if ($first) {
        $cmd = "cat $file > $of";
        $first = 0;
    } else {
        $cmd = "cat $file >> $of";
    }
    system($cmd);
    
    # check for errors (see perldoc -f system)
    if ($? == -1) {
        $logger->error("failed to execute: $!");
        exit(-1);
    }
    elsif ($? & 127) {
        my $err = sprintf("child died with signal %d, %s coredump\n", ($? & 127), ($? & 128) ? 'with':'without');
        $logger->error($err);
        exit(1);
    }
    else {
        my $exitval = $? >> 8;
        if ($exitval != 0) {
            $logger->error("child exited with value $exitval");
            exit($exitval);
        }
    }
}

if ($options->{'delete_worker_files'}) {
    my $nwf = scalar(@paths);
    $logger->info("deleting $nwf dnadist_worker output file(s)");
    foreach my $file (@paths) {
        my $count = unlink $file;
        die "error unlinking $file" if ($count != 1);
    }
}

if (($options->{'delete_tmp_files'}) && (!$options->{'worker_output_is_sorted'})) {
    my $ns = scalar(@$sorted_files);
    $logger->info("deleting $ns sorted temporary file(s)");
    foreach my $file (@$sorted_files) {
        my $count = unlink $file;
        die "error unlinking $file" if ($count != 1);
    }
}

exit(0);

## subroutines
sub check_parameters {
    my $options = shift;
    
    ## make sure required parameters were passed
    my @required = qw(worker_output_dir output_file);
    for my $option ( @required ) {
        unless  ( defined $options->{$option} ) {
            die "--$option is a required option";
        }
    }

    ## additional parameter checking
    my $owd = $options->{'worker_output_dir'};
    if (!-e $owd) {
        die "--worker_output_dir=$owd does not exist";
    }

    my $of =  $options->{'output_file'};
    if (-e $of) {
        die "--output_file=$of already exists";
    }
    
    ## defaults
    $options->{'tmp_dir'} = $DEFAULT_TMP_DIR if (!defined($options->{'tmp_dir'}));
    $options->{'worker_output_regex'} = $DEFAULT_WORKER_OUTPUT_REGEX if (!defined($options->{'worker_output_regex'}));
}

# simple in-memory sort of a dnadist_worker output file
sub sort_file {
    my($infile, $outfile) = @_;
    my $data = [];

    # read input
    my $ifh = FileHandle->new();
    $ifh->open($infile, 'r') || die "unable to read from $infile";
    my $lnum = 0;
    while (my $line = <$ifh>) {
        chomp($line);
        ++$lnum;
        if ($line =~ /^(\d+)\s+(\S+)/) {
            push(@$data, [$1, $2]);
        } else {
            die "unable to parse line $lnum of $infile: $line";
        }
    }
    $ifh->close();

    # sort and output
    my $ofh = FileHandle->new();
    $ofh->open(">$outfile") || die "unable to write to $outfile";
    foreach my $d (sort { $a->[0] <=> $b->[0] } @$data) {
        $ofh->print($d->[0] . " " . $d->[1] . "\n");
    }
    $ofh->close();
}
