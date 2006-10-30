#!/usr/local/bin/perl

## this section is just to help pull the load off the filer and
##  will probably be removed once we have a better solution.
##  for testing, you'll have to set 
##  PERL_MOD_DIR = /usr/local/devel/ANNOTATION/ard/current/lib/site_perl/5.8.5
BEGIN{foreach (@INC) {s/\/usr\/local\/packages/\/local\/platform/}};
use lib (@INC,$ENV{"PERL_MOD_DIR"});
no lib "$ENV{PERL_MOD_DIR}/i686-linux";
no lib ".";

=head1 NAME

template.pl - you should put a one-line description of your script here

=head1 SYNOPSIS

USAGE: template.pl 
            --some_argument=/path/to/some_file.fsa 
            --another_argument=/path/to/somedir
          [ --optional_argument=/path/to/somefile.list 
            --optional_argument2=1000
          ]

=head1 OPTIONS

B<--some_argument,-i>
    here you'll put a longer description of this argument that can span multiple lines. in 
    this example, the script has a long option and a short '-i' alternative usage.

B<--another_argument,-o>
    here's another named required argument.  please try to make your argument names all
    lower-case with underscores separating multi-word arguments.

B<--optional_argument,-s>
    optional.  you should preface any long description of optional arguments with the
    optional word as I did in this description.  you shouldn't use 'optional' in the
    actual argument name like I did in this example, it was just to indicate that
    optional arguments in the SYNOPSIS should go between the [ ] symbols.

B<--optional_argument2,-f>
    optional.  if your optional argument has a default value, you should indicate it
    somewhere in this description.   ( default = foo )

B<--debug> 
    Debug level.  Use a large number to turn on verbose debugging. 

B<--log,-l> 
    Log file

B<--help,-h>
    This help message

=head1  DESCRIPTION

put a longer overview of your script here.

=head1  INPUT

the input expectations of your script should be here.  pasting in examples of file format
expected is encouraged.

=head1  OUTPUT

the output format of your script should be here.  if your script manipulates a database,
you should document which tables and columns are affected.

=head1  CONTACT

    Your Name
    you@tigr.org

=cut

use strict;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Pod::Usage;
use Ergatis::Logger;

my %options = ();
my $results = GetOptions (\%options, 
                          'some_argument|i=s',
                          'another_argument|o=s',
                          'optional_argument|s=s',
                          'optional_argument2|f=s',
                          'log|l=s',
                          'debug=s',
                          'help|h') || pod2usage();

my $logfile = $options{'log'} || Ergatis::Logger::get_default_logfilename();
my $logger = new Ergatis::Logger('LOG_FILE'=>$logfile,
                                  'LOG_LEVEL'=>$options{'debug'});
$logger = $logger->get_logger();

## display documentation
if( $options{'help'} ){
    pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDERR} );
}

## make sure everything passed was peachy
&check_parameters(\%options);

##
## CODE HERE
##

exit(0);

sub check_parameters {
    my $options = shift;
    
    ## make sure required arguments were passed
    unless ( $options{some_argument} && $options{another_argument} ) {
        $logger->logdie("You must pass both --some_argument and --another_argument");
    }
    
    ##
    ## you can do other things here, such as checking that files exist, etc.
    ##
    
    ## handle some defaults
    $options{optional_argument2}   = 'foo'  unless ($options{optional_argument2});
}
