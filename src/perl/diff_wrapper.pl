#!/usr/local/bin/perl
BEGIN{foreach (@INC) {s/\/usr\/local\/packages/\/local\/platform/}};
use lib (@INC,$ENV{"PERL_MOD_DIR"});
no lib "$ENV{PERL_MOD_DIR}/i686-linux";
no lib ".";

=head1  NAME 

dummy.pl - do nothing

=head1 SYNOPSIS

USAGE:  dummy.pl --expected file.expected --result file.result --debug_level debug_level --log log_file

=head1 OPTIONS

=item *

B<--debug,-d> Debug level.  Use a large number to turn on verbose debugging. 

=item *

B<--expected> File containing the expect query results

=item *

B<--result> File containing the most recent query results

=item *

B<--log,-l> Log file

=item *

B<--help,-h> This help message

=back

=head1   DESCRIPTION

=cut


use strict;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use Pod::Usage;
BEGIN {
use Workflow::Logger;
}

my %options = ();
my $results = GetOptions (\%options,
			  'expected=s',
			  'result=s',
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

print `diff $options{'expected'} $options{'result'}`;


sub check_parameters{
    my ($options) = @_;

    my $expected;
    my $result;

    if (( exists $options{'expected'}) && (defined($options{'expected'}))){
	# nothing
    }
    else {
	pod2usage({-exitval => 2,  -message => "--expected option missing", -verbose => 1, -output => \*STDERR});    
	$logger->logdie("expected was not defined");
    }

    if (( exists $options{'result'}) && (defined($options{'result'}))){
	# nothing
    }
    else {
	pod2usage({-exitval => 2,  -message => "--result option missing", -verbose => 1, -output => \*STDERR});    
	$logger->logdie("result was not defined");
    }
}
