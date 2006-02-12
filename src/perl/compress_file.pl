#!/usr/local/bin/perl
use lib (@INC,$ENV{"PERL_MOD_DIR"});
no lib "$ENV{PERL_MOD_DIR}/i686-linux";
no lib ".";

=head1  NAME 

compress_file.pl - wrapper to gzip an output file

=head1 SYNOPSIS

USAGE: augustus2bsml.pl 
        --file=/path/to/some.file 
      [ --compress=1
        --list_file=/path/to/some.list
        --log=/path/to/some.log
        --debug=4
      ]

=head1 OPTIONS

B<--file,-f> 
    The path to any file you wish to compress.

B<--compress,-c>
    Optional.  If 0, the compression will not be performed.

B<--compress,-c>
    Optional.  Writes a list file of the gzipped file created.

B<--debug,-d> 
    Debug level.  Use a large number to turn on verbose debugging. 

B<--log,-l> 
    Log file

B<--help,-h> 
    This help message

=head1   DESCRIPTION

This script is used as a gzip wrapper to compress a file.

=head1 INPUT

The input is any file you wish to compress, usually a text file such
as alignment data.  This file is passed with the --file option.

The --compress option exists so that this script can be incorporated
into workflow.  Since Workflow is not a logic-based system, it cannot
call gzip on output *sometimes*.  So we call this script all the time,
and just pass --compress 1 or 0 depending on whether we actually want
the compression done.

The --list_file options allows us to create a list file containing the
full path to the list file created.  This is useful because a workflow
may create a list file for the raw output, but if we zip it, we want to
overwrite that single-entry list file with the gzip'ed path before it
is merged.

=head1 OUTPUT

The output file will be named the same but with the .gz extension

=head1 CONTACT

    Joshua Orvis
    jorvis@tigr.org

=cut


use strict;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use Pod::Usage;
BEGIN {
use Workflow::Logger;
}

my %options = ();
my $results = GetOptions (\%options,
			  'file|f=s',
			  'compress|c=s',
              'list_file|s=s',
			  'log|l=s',
			  'debug=s',
			  'help|h') || pod2usage();

my $logfile = $options{'log'} || Workflow::Logger::get_default_logfilename();
my $logger = new Workflow::Logger('LOG_FILE'=>$logfile,
				  'LOG_LEVEL'=>$options{'debug'});
$logger = Workflow::Logger::get_logger();

# display documentation
if( $options{'help'} ){
    pod2usage( {-exitval=>0, -verbose => 2, -output => \*STDERR} );
}

&check_parameters(\%options);

my $return_status = 0;

## we're only going to compress if the user requested it
if ($options{compress}) {
    $logger->debug("preparing to compress file $options{file}") if ($logger->is_debug);
    system("gzip -f $options{file}");
    
    ## was it successful?
    $return_status = $? >> 8;
    $logger->debug("return status was $return_status") if ($logger->is_debug);
} else {
    $logger->debug("skipping compression of file $options{file} because --compress=$options{compress}") if ($logger->is_debug);
}

if ($options{list_file}) {
    $logger->debug("creating list file $options{list_file}") if ($logger->is_debug);
    open(my $lfh, ">$options{list_file}") || $logger->logdie("can't create list file $options{list_file}");
    
    print $lfh "$options{file}";
    print $lfh ".gz" if $options{compress};
    print $lfh "\n";
} else {
    $logger->debug("skipping list file creation because --list_file=$options{list_file}") if ($logger->is_debug);
}

exit($return_status);

sub check_parameters{
    my ($options) = @_;
    
    if ($options{'file'} eq ""){
	    pod2usage({-exitval => 2,  -message => "--file option missing", -verbose => 1, -output => \*STDERR});    
    }
    
    ## make sure the file exists
    if (! -e $options{'file'}){
	    $logger->logdie("output file $options{file} doesn't exist");
    }
    
    ## handle some defaults
    $options{compress}  = 1 if (! defined $options{compress});
    $options{list_file} = 0 if (! defined $options{list_file});
}
