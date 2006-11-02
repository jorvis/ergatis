#!/usr/local/bin/perl
BEGIN{foreach (@INC) {s/\/usr\/local\/packages/\/local\/platform/}};
use lib (@INC,$ENV{"PERL_MOD_DIR"});
no lib "$ENV{PERL_MOD_DIR}/i686-linux";
no lib ".";

=head1  NAME 

trf_wrapper.pl - wraps tandem repeat finder (trf) to allow proper return values
and output location options.

=head1 SYNOPSIS

USAGE: trf2bsml.pl 
            --input=/path/to/somefile.dat 
          [ --output_dir=/path/to/output.bsml ]
            --match=2
            --mismatch=7
            --delta=7
            --pm=80
            --pi=10
            --minscore=50
            --maxperiod=500
          [ --other_opts='-d' ]

=head1 OPTIONS

B<--input,-i> 
    input .dat file from a RepeatMasker search.

B<--output_dir,-o> 
    [optional] Directory where you want the output written.  The script will chdir into
    that directory, run trf, and then chdir back.

B<--match,-m>
    matching weight.

B<--mismatch,-x>
    mismatching penalty.

B<--delta,-e>
    indel penalty.

B<--pm,-p>
    match probability (whole number)

B<--pi,-r>
    indel probability (whole number)

B<--minscore,-n>
    minimum alignment score to report.

B<--maxperiod,-a>
    maximum period size to report.

B<--other_opts,-t>
    This should be put in quotes.  Acceptable options are:
    
    -m masked sequence file
    -f flanking sequence
    -d data file

B<--debug,-d> 
    Debug level.  Use a large number to turn on verbose debugging. 

B<--log,-l> 
    Log file

B<--help,-h> 
    This help message

=head1   DESCRIPTION

Tandem Repeat Finder (trf) has some interesting quirks that makes it difficult
to run directly using workflow.  It does not return 0 upon success, it returns the
number of sequences processed, and it has no options to control the output file
names or locations.

This wrapper allows you to pass the directory where you want the out files to go,
and checks the return value of trf, logs appropriately, and returns 0 upon success.

=head1 INPUT

Each of the options available to trf can be used by passing them to this script.  See
the OPTIONS section for descriptions of these.  The input FASTA file can have one
or many sequences within it.

=head1 OUTPUT

The output files will be written in the location defined using the --output_dir option.
If not passed, /tmp will be used.  

=head1 CONTACT

    Joshua Orvis
    jorvis@tigr.org

=cut

use strict;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use Pod::Usage;
BEGIN {
use Ergatis::Logger;
}

my %options = ();
my $results = GetOptions (\%options,
              'input|i=s',
              'output_dir|o=s',
			  'match|m=s',
              'mismatch|x=s',
              'delta|e=s',
              'pm|p=s',
              'pi|r=s',
              'minscore|n=s',
              'maxperiod|a=s',
              'other_opts|t=s',
              'debug|d=s',
              'log|l=s',
			  'help|h') || pod2usage();

my $logfile = $options{'log'} || Ergatis::Logger::get_default_logfilename();
my $logger = new Ergatis::Logger('LOG_FILE'=>$logfile,
				  'LOG_LEVEL'=>$options{'debug'});
$logger = $logger->get_logger();

# display documentation
if( $options{'help'} ){
    pod2usage( {-exitval=>0, -verbose => 2, -output => \*STDOUT} );
}

## make sure all passed options are peachy
&check_parameters(\%options);

## get the current working directory
my $starting_dir = `pwd`;

## cd to the run area
chdir($options{output_dir}) || $logger->logdie("can't cd to temp space $options{output_dir}");

## run the command
my $cmd = "/usr/local/bin/trf $options{input} $options{match} $options{mismatch} $options{delta} $options{pm} $options{pi} $options{minscore} $options{maxperiod} $options{other_opts}";
system($cmd);

## from the author: trf returns the number of sequences processed.  failures are reported
##  with a negative return code.
if ($? < 1) {
    $logger->logdie("trf failed returning $?");
}

## go back to the original directory
chdir($starting_dir);

exit(0);

sub check_parameters {
    
    ## make sure input file exists
    if (! -e $options{'input'}) { $logger->logdie("input file $options{'input'} does not exist") }
   
    ## handle defaults
    $options{output_dir}   = '/tmp/' unless ($options{output_dir});
    $options{other_opts} = ''      unless ($options{other_opts});
    
    return 1;
}
