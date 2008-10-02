#!/usr/bin/perl

=head1  NAME 

dnadist_combiner_wrapper.pl - Wrapper for dnadist_combiner executable in a parallel version of dnadist.

=head1 SYNOPSIS

dnadist_combiner_wrapper.pl
         --dnadist_combiner_path=/path/to/dnadist_combiner
         --control_file=/path/to/dnadist_control_file
         --sorted_worker_output=/path/to/partial_output.dnadist
         --output_file=/path/to/output.txt
        [--log=/path/to/some.log
         --debug=4 ]

=head1 OPTIONS

B<--dnadist_combiner_path,-w>
    path to the dnadist_combiner executable.

B<--control_file,-c>
    path to the dnadist control file.  see the dnadist documentation for details (this
    file is essentially a set of keyboard commands that gets piped into dnadist's 
    interactive menu interface)

B<--sorted_worker_output,-s>
    path to the sorted, concatenated output from all of the distributed dnadist_worker
    processes.  dnadist_combiner will fail if this output is not complete or does not
    match the input file passed in the --control_file

B<--output_file,-o>
    path to the file in which to store the final dnadist output.

B<--log,-l> 
    optional.  path to a log file the script should create.  will be overwritten if
    already exists.
    
B<--debug,-d> 
    optional.  the debug level for the logger (an integer)

B<--help,-h> 
    This help message/documentation.

=head1 DESCRIPTION

Wrapper for dnadist_combiner executable in a parallel version of dnadist.

=head1 INPUT

The input to this wrapper consists of the input to the original dnadist
program, specified by the --control_file, plus the tab-delimited data
produced by the dnadist_worker processes.  Each line of this output will
have two columns separated by a space, e.g.

1 0.244838
2 0.056921
3 0.647833
4 0.234556

This file contains the sorted concatenated output from all of the distributed
dnadist_worker jobs.

=head1 OUTPUT

The output of the wrapper is the output of dnadist, namely a distance matrix
whose format is specified by the --control_file.

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
            "dnadist_combiner_path|w=s",
            "control_file|c=s",
            "sorted_worker_output|s=s",
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

# construct dnadist_combiner command line
my $cmd = "cat " . $options->{'control_file'} . " | ";
$cmd .= join(' ', map {$options->{$_}} qw(dnadist_combiner_path sorted_worker_output output_file));
$logger->debug("running " . $cmd);

# run dnadist_combiner
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
    my @required = qw(dnadist_combiner_path control_file sorted_worker_output output_file);
    for my $option ( @required ) {
        unless  ( defined $options->{$option} ) {
            die "--$option is a required option";
        }
    }

    ## additional parameter checking
    my $wp =  $options->{'dnadist_combiner_path'};
    if (!-e $wp || !-x $wp) {
        die "$wp does not exist or is not executable";
    }

    my $cf =  $options->{'control_file'};
    if (!-e $cf) {
        die "$cf does not exist";
    }

    my $wf = $options->{'sorted_worker_output'};
    if (!-e $wf) {
        die "--sorted_worker_output=$wf does not exist";
    }

    my $of =  $options->{'output_file'};
    if (defined($of) && (-e $of)) {
        die "--output_file=$of already exists";
    }
    
    ## defaults
}
