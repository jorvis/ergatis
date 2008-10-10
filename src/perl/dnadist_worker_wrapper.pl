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
        [--instance_limit=2
         --log=/path/to/some.log
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

B<--instance_limit,-i>
    optional.  run no more than this number of dnadist_worker_wrapper jobs on any
    given machine (this is an ad-hoc architecture-dependent workaround; use at
    your peril.)

B<--log,-l> 
    optional.  path to a log file the script should create.  will be overwritten if
    it already exists.
    
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
            "instance_limit|i=i",
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
&enforce_instance_limit($logger, $options->{'instance_limit'});

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

    if (defined($options->{'instance_limit'})) {
        if ($options->{'instance_limit'} < 1) {
            $options->{'instance_limit'} = 1;
        }
    }
    
    ## defaults
}

# this is an ad-hoc workaround for problems introduced when too many dnadist_worker
# processes are active on the same machine.
sub enforce_instance_limit {
    my($logger, $limit) = @_;

    if (defined($limit)) {
        # always be extra-careful before entering an infinite loop...
        $limit = 1 if (($limit < 1) || ($limit !~ /^\d+$/));
        my $self_regex = $0;
        $self_regex =~ s/\//\\\//g;

        while (1) {
            my $ps = `ps -f`;
            my $procs = [];
            my $ppid = undef;
            my $my_ppid = undef;
            
            $logger->debug("$$: checking 'ps' output");
            foreach my $proc (split(/\n/, $ps)) {
                my($uid, $pid, $ppid, $cmd) = ($proc =~ /^(\d+)\s+(\d+)\s+(\d+)\s+\d+\s+[\d\:]+\s+\S+\s+\S+\s+(\S.*)$/);
                if ($cmd =~ /$self_regex/) {
                    $logger->debug("$$: uid=$uid pid=$pid ppid=$ppid cmd=$cmd");
                    push(@$procs, [$uid, $pid, $ppid, $cmd]);
                    $my_ppid = $ppid if ($pid == $$);
                }
            }
    
            my @siblings = grep { $_->[2] == $my_ppid } @$procs;
            my @sorted = sort { $a->[1] <=> $b->[1] } @siblings;

            my $rank = undef;
            my $ns = scalar(@sorted);
            for (my $i = 0;$i < $ns;++$i) {
                if ($sorted[$i]->[1] == $$) {
                    $rank = $i + 1;
                    last;
                }
            }

            $logger->debug("$$: my_ppid=$my_ppid, rank=$rank, limit=$limit");
            last if ($rank <= $limit);
            $logger->debug("$$: sleeping");
            sleep(30);
        }
        $logger->debug("$$: running");
    }
}
