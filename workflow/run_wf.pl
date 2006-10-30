#!/local/perl/bin/perl
use lib (@INC,$ENV{"PERL_MOD_DIR"});
no lib "$ENV{PERL_MOD_DIR}/i686-linux";
no lib ".";

=head1  NAME 

run_wf.pl - do nothing

=head1 SYNOPSIS

USAGE:  run_wf.pl 
    --instance instance file 
    --template xml template file 
    --ini ini file 
  [ --skiprun
    --nodistib
	--workflow_root root of workflow install
    --debug debug_level
    --log log_file
  ]

=head1 OPTIONS

=item *

B<--debug,-d> Debug level.  Use a large number to turn on verbose debugging. 

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
use Ergatis::Logger;
use Ergatis::Run;

my %options = ();
my $results = GetOptions (\%options, 
			  'instance|i=s',
			  'ini|c=s',
			  'template|t=s',
			  'skiprun=s',
			  'nodistrib=s',
  			  'workflow_root:s',
			  'log|l=s',
			  'debug=s',
			  'help|h') || pod2usage();

## friendly permissions
umask(0000);

my $logfile = $options{'log'} || Ergatis::Logger::get_default_logfilename();
my $logger = new Ergatis::Logger('LOG_FILE'=>$logfile,
				  'LOG_LEVEL'=>$options{'debug'});
$logger = Ergatis::Logger::get_logger();

# display documentation
if( $options{'help'} ){
    pod2usage( {-exitval=>0, -verbose => 2, -output => \*STDERR} );
}

&check_parameters(\%options);

my $wfexec;
if ($options{'workflow_root'}) {
	$wfexec = new Ergatis::Run( 'nodistrib'     => $options{'nodistrib'},
								 'WORKFLOW_ROOT' => $options{'workflow_root'}
							   );
} else {
	$wfexec = new Ergatis::Run('nodistrib'=>$options{'nodistrib'});
}

my $instancexmlfile= "$options{'instance'}";

if(!(-e $instancexmlfile)){
    my $createlogfile = "$instancexmlfile.create.log";
    my $createoutfile = "$instancexmlfile.create.out";

    my $createstatus = $wfexec->CreateWorkflow($instancexmlfile, $options{'ini'}, $options{'template'}, $createlogfile, $createoutfile);
    if($createstatus == 0){
	$logger->debug("Created instance file $instancexmlfile") if($logger->is_debug());
    }
    else{
	$logger->logdie("Unable to create instance file $instancexmlfile");
	exit ($createstatus);
    }
}
if(!($options{'skiprun'})){
    my $runlogfile = "$instancexmlfile.run.log";
    my $runoutfile = "$instancexmlfile.run.out";
    
    my $runstatus = $wfexec->RunWorkflow($instancexmlfile,$runlogfile,$runoutfile);
    if($runstatus == 0){
	$logger->debug("Completed execution of instance file $instancexmlfile") if($logger->is_debug());
	print "Running instance file $instancexmlfile\n";
	exit(0);
    }
    else{
	$logger->logdie("Unable to run instance file $instancexmlfile. Return code $runstatus");
	exit ($runstatus);
    }
}
else{
    $logger->get_logger()->info("Skipping execution of workflow $instancexmlfile");
	print "Skipping execution of workflow $instancexmlfile\n";
}

sub check_parameters{
    my ($options) = @_;
    
    if(0){
	pod2usage({-exitval => 2,  -message => "error message", -verbose => 1, -output => \*STDERR});    
    }
}
