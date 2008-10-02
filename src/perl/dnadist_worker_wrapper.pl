#!/usr/bin/perl

=head1  NAME 

dnadist_worker_wrapper.pl - Wrapper for dnadist_worker executable in a parallel version of dnadist.

=head1 SYNOPSIS

dnadist_worker_wrapper.pl
         --dnadist_worker_path=/path/to/dnadist_worker
         --control_file=/path/to/dnadist_control_file
         --start_cell=10
         --end_cell=20
         --output_file=/path/to/partial_output.dnadist
        [--log=/path/to/some.log
         --debug=4 ]

=head1 OPTIONS

B<--dnadist_worker_path,-w>
    path to the dnadist_worker executable.

B<--control_file,-c>
    path to the dnadist control file.  see the dnadist documentation for details (this
    file is essentially a set of keyboard commands that gets piped into dnadist's 
    interactive menu interface)

B<--start_cell,-s> 
    first matrix cell (inclusive) for the worker process to calculate

B<--end_cell,-e> 
    last matrix cell (inclusive) for the worker process to calculate

B<--output_file,-o>
    path to the file in which to store the (partial) dnadist output.

B<--log,-l> 
    optional.  path to a log file the script should create.  will be overwritten if
    already exists.
    
B<--debug,-d> 
    optional.  the debug level for the logger (an integer)

B<--help,-h> 
    This help message/documentation.

=head1 DESCRIPTION

Wrapper for dnadist_worker executable in a parallel version of dnadist.

=head1 INPUT

The input to this wrapper is the same as the input to the original dnadist
program, and is specified by the --control_file.

=head1 OUTPUT

The (partial) dnadist output will be written to the file specified by --output_file.
Each line of the output will have two columns separated by a space, e.g.

1 0.244838
100 0.244838
2 0.056921
200 0.056921

The first value on each line is the 0-based offset of a cell in the output
matrix and the second value is the distance value computed for that cell.
Out-of-order offsets appear in the file because unless dnadist is configured
to output a lower-triangular matrix each distance value will appear in two 
distinct (but symmetric) locations.

=head1 CONTACT

    Jonathan Crabtree
    jonathancrabtree@gmail.com

=cut

use strict;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Pod::Usage;
use Ergatis::Logger;

## input
my $options = {};
&GetOptions($options,
            "dnadist_worker_path|w=s",
            "control_file|c=s",
            "start_cell|s=i",
            "end_cell|e=i",
            "output_file|o=s",
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

# construct dnadist_worker command line
my $cmd = "cat " . $options->{'control_file'} . " | ";
$cmd .= join(' ', map {$options->{$_}} qw(dnadist_worker_path start_cell end_cell output_file));
$logger->debug("running " . $cmd);

# run dnadist_worker
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

exit(0);

## subroutines
sub check_parameters {
    my $options = shift;
    
    ## make sure required parameters were passed
    my @required = qw(dnadist_worker_path control_file start_cell end_cell output_file);
    for my $option ( @required ) {
        unless  ( defined $options->{$option} ) {
            die "--$option is a required option";
        }
    }

    ## additional parameter checking
    my $wp =  $options->{'dnadist_worker_path'};
    if (!-e $wp || !-x $wp) {
        die "$wp does not exist or is not executable";
    }

    my $cf =  $options->{'control_file'};
    if (!-e $cf) {
        die "$cf does not exist";
    }

    my $of =  $options->{'output_file'};
    if (-e $of) {
        die "--output_file=$of already exists";
    }

    if ($options->{'start_cell'} > $options->{'end_cell'}) {
        die "--start_cell must be <= --end_cell";
    }
    
    ## defaults
}
