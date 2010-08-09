#!/usr/bin/perl

=head1 NAME

run_genemark-es.pl - Executes gm_es.pl from the GeneMark-ES package but with magic around it
to make up for some 'problems' with the original program.

=head1 SYNOPSIS

USAGE: run_cmsearch.pl 
              --input_file=/path/hmmpfam/result.raw
              --tmp_dir=/path/to/tmpdir/
              --output_file=/path/to/infernal.raw
              --other_opts='foo bar'
              --executable=/path/to/gm_es.pl
              --log=/path/to/some/file.log

=head1 OPTIONS

B<--input_file,-i>
    An input, single-sequence FASTA file

B<--tmp_dir,-t>
    Directory to write temporary files.  Must exist already, and will not be cleaned up.

B<--output_file,-o>
    The output file for results to go into.

B<--other_opts,-e>
    Other options to be passed into the GeneMark-ES program.  

B<--executable,-x>
    The path to the executable that will be run.  

B<--debug> 
    Debug level.  Use a large number to turn on verbose debugging. 

B<--log,-l> 
    Log file

B<--help,-h>
    This help message

=head1  DESCRIPTION

The genemark-es script dumps temporary files and directories all over your current
working directory.  This wrapper moves into a 'safe' directory, executes the script
from there, then cleans up after itself.

=head1  INPUT


=head1  OUTPUT


=head1  CONTACT

    Joshua Orvis
    jorvis@users.sf.net

=cut

use strict;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Pod::Usage;
use Ergatis::Logger;
use File::Copy;

my %options = ();
my $results = GetOptions (\%options, 
                          'input_file|i=s',
                          'tmp_dir|t=s',
                          'output_file|o=s',
                          'other_opts|e=s',
                          'executable|x=s',
                          'log|l=s',
                          'debug=s',
                          'help|h') || &_pod();

## display documentation
&_pod if( $options{'help'} );

my $logfile = $options{'log'} || Ergatis::Logger::get_default_logfilename();
my $logger = new Ergatis::Logger('LOG_FILE'=>$logfile,
				  'LOG_LEVEL'=>$options{'debug'});
$logger = $logger->get_logger();

## make sure everything passed was peachy
&check_parameters(\%options);

my $command = "$options{executable} ";

if ( $options{other_opts} ) {
    $command .= "$options{other_opts} ";
}

$command .= $options{input_file};

chdir( $options{tmp_dir} );

if ( system($command) ) {
    _log("failed to run command ($command): returned $?");
}

## copy the output file
copy('genemark_hmm.gtf', $options{output_file}) || die "failed to copy output file to $options{output_file}: $!";

## now clean up
for ( qw(data dna.fa.X info mod pred run) ) {
    cleanup($_);
}

exit(0);

sub check_parameters {
    my $options = shift;
    
    ## make sure required arguments were passed
    my @required = qw( input_file tmp_dir output_file executable );
    for my $option ( @required ) {
        unless  ( defined $$options{$option} ) {
            die "--$option is a required option";
        }
    }
    
    ##
    ## you can do other things here, such as checking that files exist, etc.
    ##
    
    ## handle some defaults
    #$options{optional_argument2}   = 'foo'  unless ($options{optional_argument2});
}

sub cleanup {
    my $dir = shift;
	local *DIR;

	opendir DIR, $dir or die "failed to read directory $dir: $!";
    
    for (readdir DIR) {
        next if /^\.{1,2}$/;
        my $path = "$dir/$_";
        unlink $path if -f $path;
        cleanup($path) if -d $path;
	}
	closedir DIR;
	rmdir $dir or print "error - $!";
}

























