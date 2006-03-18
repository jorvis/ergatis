#!/usr/local/bin/perl
BEGIN{foreach (@INC) {s/\/usr\/local\/packages/\/local\/platform/}};
use lib (@INC,$ENV{"PERL_MOD_DIR"});
no lib "$ENV{PERL_MOD_DIR}/i686-linux";
no lib ".";
#------------------------------------------------------------------------------
#
# author:   Jay Sundaram sundaram@tigr.org
#
# date:     2005-12-14
#
# cvs:      ergatis/src/perl/gatekeeper.pl
#
# $Id$
#
#
#------------------------------------------------------------------------------
=head1 NAME

gatekeeper.pl - Creates/removes database lock file from project repository root

=head1 SYNOPSIS

USAGE:  gatekeeper.pl -D database [-U username] -a action -c component [-d debug_level] [-h] [-l log4perl] [-m] -r repository -p pipeline

=head1 OPTIONS

=over 8

=item B<--database,-D>

    Name of the database to be locked

=item B<--action,-a>

    Action to be taken either 'create' or 'remove'

=item B<--debug_level,-d>

    Optional: Coati::Logger log4perl logging level.  Default is 0

=item B<--help,-h>

    Print this help

=item B<--log4perl,-l>

    Optional - log4perl log file.  Default is /tmp/gatekeeper.pl.log

=item B<--man,-m>

    Display the pod2usage page for this utility

=item B<--pipeline,-p>

    The fullpath to the current pipeline.xml.instance file

=item B<--repository,-r>

    The project's repository root

=item B<--username,-U>

    Optional - Username of person to be notified via email. Default is value returned from whoami

=item B<--component,-c>

    Name of the database manipulating workflow component e.g. initdb or legacy2bsml

=back

=head1 DESCRIPTION

    gatekeeper.pl - Creates/removes database lock file from project repository root

    Assumptions:
    1. User has appropriate permissions (to execute script, access chado database, write to output directory).
    2. All software has been properly installed, all required libraries are accessible.

    Sample usage:
   ./gatekeeper.pl -D chado_test -a create -l gatekeeper.pl.log -U sundaram -r /usr/local/annotation/CHADO_TEST -c bsml2chado -p /usr/local/annotation/CHADO_TEST/Workflow/pipeline/376/pipeline.xml.instance

=head1 CONTACT

    Jay Sundaram
    sundaram@tigr.org

=cut

use Mail::Mailer;
use strict;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use Pod::Usage;
use Data::Dumper;
use XML::Twig;

BEGIN {
use Workflow::Logger;
}


$|=1;

#-------------------------------------------------------------
# Parse command line options
#-------------------------------------------------------------

my ($debug_level, $help, $log4perl, $man, $pipeline, $username, $repository, $component, $action, $database);


my $results = GetOptions (
			  'action|a=s'       => \$action,
			  'log4perl|l=s'     => \$log4perl,
			  'debug_level|d=s'  => \$debug_level, 
			  'help|h'           => \$help,
			  'man|m'            => \$man,
			  'pipeline|p=s'     => \$pipeline,
			  'username|U=s'     => \$username,
			  'repository|r=s'   => \$repository,
			  'database|D=s'     => \$database,
			  'component|c=s'    => \$component
			  );


&pod2usage({-exitval => 1, -verbose => 2, -output => \*STDOUT}) if ($man);
&pod2usage({-exitval => 1, -verbose => 1, -output => \*STDOUT}) if ($help);


$debug_level = 5;

#
# initialize the logger
#
if (!defined($log4perl)){
    $log4perl = "/tmp/gatekeeper.pl.log";
}

my $mylogger = new Workflow::Logger('LOG_FILE'=>$log4perl,
				    'LOG_LEVEL'=>$debug_level);

my $logger = Workflow::Logger::get_logger(__PACKAGE__);

my $errorctr=0;

if (!defined($database)){
    print STDERR "database was not defined\n";
    $errorctr++;
}


if (!defined($action)){
    print STDERR "action was not defined\n";
    $errorctr++;
}

if (!defined($component)){
    print STDERR "component was not defined\n";
    $errorctr++;
}

if (!defined($repository)){
    print STDERR "repository was not defined\n";
    $errorctr++;
}

if (!defined($pipeline)){
    print STDERR "pipeline was not defined\n";
    $errorctr++;
}

if ($errorctr > 0 ) {
    &print_usage();
}



#
# username must be specified
#
if (!defined($username)){
    $username = `whoami`;
    chomp $username;
    $logger->warn("username was set as '$username'");
}

#
# database must be specified
#
if (!defined($database)){
    $logger->logdie("database was not defined");
}

#
# action must be specified
#
&verify_action($action);

#
# repository must be specified
#
&verify_repository_root($repository);


my $lockfile = $repository . "/database_lock";


if (-e $lockfile){

    #
    # The database lock file exists
    #
    if ($action eq 'create'){
	

	if (&pipeline_running($lockfile)){
	    
	    &notify_user( $lockfile,
			  $username,
			  $database,
			  $repository,
			  $pipeline,
			  $component,
			  $action );
	    
	    #
	    # Abort execution of the pipeline by having this script die 
	    #
	    $logger->logdie("Found database lock file '$lockfile'. Please review logfile '$log4perl'");
	}
	else {

	    unlink $lockfile;

	    &create_lock_file($lockfile,
			      $username,
			      $database,
			      $repository,
			      $pipeline,
			      $component,
			      $action );

	}
	    
    }
    elsif ($action eq 'remove') {
	
	if (&pipeline_running($lockfile, $pipeline, $action)){


	    $logger->logdie("Cannot remove database lock file '$lockfile' because a database manipulating workflow is running.  action '$action' username '$username' database '$database' repository '$repository' component '$component' pipeline '$pipeline'");

	}
	else {
	    $logger->warn("Removing the following database lock file '$lockfile'.");
	    $logger->warn("username [$username] database [$database] repository [$repository] pipeline [$pipeline] component [$component] action [$action]");
	    
	    unlink $lockfile;
	}
    }
    else {
	$logger->logdie("Unrecognized action '$action'");
    }
}
else {
    #
    # The database lock file does not exist
    #
    if ($action eq 'create'){
	
	&create_lock_file($lockfile,
			  $username,
			  $database,
			  $repository,
			  $pipeline,
			  $component,
			  $action );
    }
    elsif ($action eq 'remove'){

	$logger->logdie("Cannot remove database lock file '$lockfile' because file does not exist.  action '$action' username '$username' database '$database' repository '$repository' component '$component' pipeline '$pipeline'");


    }
    else {
	$logger->logdie("Unrecognized action '$action'");
    }
}

#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#
#                                  END OF MAIN  -- SUBROUTINES FOLLOW
#
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#---------------------------------------------------------------------
# create_lock_file()
#
# input:
#
# output:
#
# return:
#
# comment:
#
#---------------------------------------------------------------------
sub create_lock_file {
    
    my ($lockfile, $username, $database, $repository, $pipeline, $component, $action) = @_;
    
    $logger->warn("Attempting to create database lock file '$lockfile'");

    open (OUTFILE, ">$lockfile") || $logger->logdie("Could not open database lock file '$lockfile': $!");
    
    my $date = `date`;

    chomp $date;

    my $content = "Creating a database lock file\n".
    "date [$date]\n".
    "database [$database]\n".
    "username [$username]\n".
    "component [$component]\n".
    "repository [$repository]\n".
    "action [$action]\n".
    "pipeline [$pipeline]\n";
    
    
    print OUTFILE $content;
    
}



#------------------------------------------------------
# send_notification()
#
#------------------------------------------------------
sub send_notification {

    my ($username, $subject, $body) = @_;

    if ($username !~ /\@/){
	$username .= "\@tigr.org";
    }

    my $mailer = Mail::Mailer->new ('sendmail');

    $mailer->open({
	             To      => $username,
		     Subject => $subject
		 }) or $logger->logdie("Could not create and send message");
    
    print $mailer $body;
    
    $mailer->close;
    


    if ($logger->is_debug()){
	$logger->debug("Notification sent to $username");
    }

}


#-----------------------------------------------------------------
# verify_action()
#
# input:
#
# output:
#
# return:
#
# comment:
#
#
#-----------------------------------------------------------------
sub verify_action {

    my ($action) = @_;

    if (!defined($action)){
	$logger->logdie("action was not defined");
    }
    else {
	if (($action eq 'create') || ($action eq 'remove')){
	    if ($logger->is_debug()){
		$logger->debug("action was specified as '$action'");
	    }
	}
    }

}

#--------------------------------------------------------------
# verify_repository_root()
#
# input:
#
# output:
#
# return:
#
# comment:
#
#
#--------------------------------------------------------------
sub verify_repository_root {

    my ($repository) = @_;

    if (!defined($repository)){
	$logger->logdie("repository was not defined");
    }
    else {
	if (!-e $repository){
	    $logger->logdie("repository '$repository' does not exist");
	}
	else{
	    if (!-d $repository){
		$logger->logdie("repository '$repository' is not a directory");
	    }
	    else {
		if ($logger->is_debug()){
		    $logger->debug("repository was specified as '$repository'");
		}
	    }
	}
    }
}


#-----------------------------------------------------------
# notify_user()
#
# input:
#
# output:
#
# return:
#
# comment:
#
#-----------------------------------------------------------
sub notify_user {

    my ($lockfile, $username, $database, $repository, $pipeline, $component, $action) = @_;

    my $subject = "[$database] locked - gatekeeper.pl";

    open (INFILE, "<$lockfile") || $logger->logdie("Could not open lock file '$lockfile': $!");

    my $body = "Found a database lock file [$lockfile] for\n".
    "database [$database]\n".
    "repository [$repository].\n\n".
    "Here are the contents:\n\n";

    while (my $line = <INFILE>){

	chomp $line;
	
	$body .= $line . "\n";
    }

    $body .= "\nYour pipeline [$pipeline] component [$component] shall now be aborted.  Please attempt restart your pipeline once the above mentioned pipeline has completed execution.";

    
    &send_notification( $username,
			$subject,
			$body );



    $logger->fatal($subject);

    $logger->fatal($body);

}

#------------------------------------------------------
# print_usage()
#
#------------------------------------------------------
sub print_usage {

    print STDERR "SAMPLE USAGE:  $0 -D database [-U username] -a action -c component [-d debug_level] [-h] [-l log4perl] [-m] -r repository -p pipeline_id\n".
    "  -D|--database            = Name of database\n".
    "  -U|--username            = Optional - Username of person to be notified by email (default is result of whoami)\n".
    "  -a|--action              = Action to be taken either 'create' or 'remove'\n".
    "  -c|--component           = Name of workflow component\n".
    "  -d|--debug_level         = Optional - Coati::Logger log4perl logging level.  Default is 0\n".
    "  -h|--help                = Optional - Display pod2usage help screen\n".
    "  -l|--log4perl            = Optional - Log4perl log file (default: /tmp/gatekeeper.pl.log)\n".
    "  -m|--man                 = Optional - Display pod2usage pages for this utility\n".
    "  -w|--repository          = The project's repository root\n".
    "  -p|--pipeline_id         = The pipeline.xml.instance\n";
    exit 1;

}




#------------------------------------------------------
# pipeline_running()
#
#------------------------------------------------------
sub pipeline_running {


    my ($file, $pipeline, $action) = @_;

    my $pipeline_file = &get_pipelinefile($file);



    ## if the pipeline.xml.instance exists, just process it
    if (-e "$pipeline_file.instance" ) {
	$pipeline_file .= '.instance';
    } 
    elsif ( -e "$pipeline_file.instance.gz" ) {
	$pipeline_file .= '.instance.gz';
    }
    
    my $ifh;
    if ($pipeline_file =~ /\.gz/) {
	open($ifh, "<:gzip", "$pipeline_file") || die "can't read $pipeline_file: $!"; 
    } else {
	open($ifh, "<$pipeline_file") || die "can't read $pipeline_file: $!";       
    }

    my $twig = new XML::Twig;
    $twig->parse($ifh);

    my $commandSetRoot = $twig->root;
    my $commandSet = $commandSetRoot->first_child('commandSet');

    next if (! $commandSet );
    
    my $state;

    if ( $commandSet->has_child('state') ) {
	$state  = $commandSet->first_child('state')->text();
    }
 

    if (($state eq 'complete') || 
	($state eq 'failed') ||
	($state eq 'error')){
	return 0;
    }
    else {
	if (($pipeline eq $pipeline_file) &&
	    ($action eq 'remove')){
	    #
	    # It may be the remove_database_lockfile step of the same
	    # pipeline.
	    return 0;
	}
	else {
	    return 1;
	}
    }


}



#------------------------------------------------------
# get_pipelinefile()
#
#------------------------------------------------------
sub get_pipelinefile {

    my ($file) = @_;

    my $pipelinefile;

    open (INFILE, "<$file") || $logger->logdie("Could not open database_lock file '$file': $!");

    while (my $line = <INFILE>){
	
	chomp $line;

	if ($line =~ /^pipeline\s+\[(\S+)\]/){
	    $pipelinefile = $1;
	    return $pipelinefile;
	}

    }
 
    $logger->logdie("Could not find pipeline information in database_lock file '$file'");

}



