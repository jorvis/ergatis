#!/usr/local/bin/perl
BEGIN{foreach (@INC) {s/\/usr\/local\/packages/\/local\/platform/}};
use lib (@INC,$ENV{"PERL_MOD_DIR"});
no lib "$ENV{PERL_MOD_DIR}/i686-linux";
no lib ".";
#------------------------------------------------------------------------------
#
# author:   Jay Sundaram sundaram@tigr.org
#
# date:     2006-04-17
#
# cvs:      ergatis/src/perl/notify_database_users.pl
#
# $Id$
#
#
#------------------------------------------------------------------------------
=head1 NAME

gatekeeper.pl - Notifies all users logged into specified database

=head1 SYNOPSIS

USAGE:  gatekeeper.pl -D database -U username -P password -a action -c component [-d debug_level] [-h] [-l log4perl] [-m] -r repository -p pipeline

=head1 OPTIONS

=over 8

=item B<--database,-D>

    Name of the database to be locked

=item B<--action,-a>

    Action to be taken either 'logout' or 'login'

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

    Login name 

=item B<--password,-P>

    Login password

=item B<--component,-c>

    Name of the database manipulating workflow component e.g. initdb or legacy2bsml

=back

=head1 DESCRIPTION

    notify_database_users.pl - Creates/removes database lock file from project repository root

    Assumptions:
    1. User has appropriate permissions (to execute script, access chado database, write to output directory).
    2. All software has been properly installed, all required libraries are accessible.

    Sample usage:
   ./notify_database_users.pl -D chado_test -a create -l notify_database_users.pl.log -U sundaram -r /usr/local/annotation/CHADO_TEST -c bsml2chado -p /usr/local/annotation/CHADO_TEST/Workflow/pipeline/376/pipeline.xml.instance

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

my ($debug_level, $help, $log4perl, $man, $pipeline, $username, $password, $repository, $component, $action, $database);


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
			  'component|c=s'    => \$component,
			  'password|P=s'     => \$password
			  );


&pod2usage({-exitval => 1, -verbose => 2, -output => \*STDOUT}) if ($man);
&pod2usage({-exitval => 1, -verbose => 1, -output => \*STDOUT}) if ($help);


my $logger = &set_logger($debug_level, $log4perl);

&check_arguments($username, $password, $database, $action, $component, $repository, $pipeline);

my $lockfile = $repository . "/notify_database_users";


if ($action eq 'logout'){

    my $logged_in_users = &get_logged_in_users($database, $username, $password);

    &notify_users_logout($logged_in_users, $lockfile, $database, $action, $component, $pipeline);
}
elsif ($action eq 'login'){

    my ($logged_out_users, $contents) = &get_logged_out_users_from_lockfile($lockfile);

    &notify_users_login($logged_out_users, $lockfile, $database, $action, $component, $pipeline, $contents);

}
else {
    $logger->logdie("Unexpected action '$action'");
}

print "End of program $0\n";
exit(0);




#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#
#                                  END OF MAIN  -- SUBROUTINES FOLLOW
#
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


#---------------------------------------------------------------------
# get_logged_in_users()
#
#---------------------------------------------------------------------
sub get_logged_in_users{

    my ($database, $username, $password) = @_;

    my $users = [];

    if (1){
	
	push(@{$users}, 'sundaram');

    }

    # Need to implement real functionality


    if (scalar(@{$users}) > 0) {
	return $users;
    }
    else {
	# no users
	exit(0);
    }

}


#---------------------------------------------------------------------
# get_logged_out_users_from_lockfile()
#
#---------------------------------------------------------------------
sub get_logged_out_users_from_lockfile {

    my ($file) = @_;

    if (!defined($file)){
	$logger->logdie("notify_database_users file '$file' does not exist");
    }
    if (-e $file){
	if (!-r $file){
	    $logger->logdie("notify_database_users file '$file' does not have read permissions");
	}
	
	#
	# Read contents, build user list.
	#
	open (INFILE, "<$file") or $logger->logdie("Could not open notify_database_users file '$file':$!");

	my $contents;

	while (my $line = <INFILE>){

	    $contents .= $line;
	    
	    chomp $line;
	    
	    if ($line =~ /^users \[(\S+)\]/){

		my $users = $1;
		
		my @userlist = split(/,/,$users);
		
		if (scalar(@userlist) > 0) {
		    
		    return (\@userlist, $contents);
		}
		else {
		    $logger->logdie("No users in notify_database_users file '$lockfile'");
		}
	    }
	}

	$logger->logdie("users was not defined.  Check notify_database_users file '$file' for [users] section.");
    }
    else {
	# No one was logged in.  No one to notify now.
	$logger->info("notify_database_users file does not exist");
	exit(0);
    }

}


#---------------------------------------------------------------------
# notify_users_login()
#
#---------------------------------------------------------------------
sub notify_users_login {
    
    my ($logged_out_users, $lockfile, $database, $action, $component, $pipeline, $content) = @_;

    my $users = &get_qualified_email_address_list($logged_out_users);

    my $subject = "[notify_database_users] Database '$database' is now available";

     $content = "Database '$database' in now available for general use.\n\nThe content of the notify_database_users file was:\n\n" . $content;

    &send_notification($users, $subject, $content);

    unlink $lockfile;
    
}


#---------------------------------------------------------------------
# notify_users_logout()
#
#---------------------------------------------------------------------
sub notify_users_logout {
    
    my ($logged_in_users, $lockfile, $database, $action, $component, $pipeline) = @_;

    my $content = &write_lockfile($logged_in_users, $lockfile, $database, $action, $component, $pipeline);

    my $users = &get_qualified_email_address_list($logged_in_users);

    my $subject = "[notify_database_users] Please logout of database '$database'";

    $content = "A database manipulating workflow is executing.  Please log out of database '$database' ASAP.\n\n" . $content;

    &send_notification($users, $subject, $content);
    
}
    
#---------------------------------------------------------------------
# write_lockfile()
#
#---------------------------------------------------------------------
sub write_lockfile {

    my ($logged_in_users, $lockfile, $database, $action, $component, $pipeline) = @_;

    $logger->warn("Attempting to create notify_database_users control file '$lockfile'");
    
    open (OUTFILE, ">$lockfile") || $logger->logdie("Could not open notify_database_users control file '$lockfile': $!");
    
    my $date = `date`;

    chomp $date;

    my $users = join(",", @{$logged_in_users});

    my $content = "Creating a notify_database_users control file\n\n".
    "date [$date]\n".
    "database [$database]\n".
    "component [$component]\n".
    "notify_database_users [$lockfile]\n".
    "pipeline [$pipeline]\n".
    "users [$users]\n";

    print OUTFILE $content;

    return $content;
}


#------------------------------------------------------
# send_notification()
#
#------------------------------------------------------
sub send_notification {

    my ($users, $subject, $body) = @_;

    my $mailer = Mail::Mailer->new ('sendmail');

    $mailer->open({
	             To      => $users,
		     Subject => $subject,
		     From    => undef
		 }) or $logger->logdie("Could not create and send message");
    

    print $mailer $body;
    
    $mailer->close;
    
}


#-----------------------------------------------------------------
# get_qualified_email_address_list()
#
#-----------------------------------------------------------------
sub get_qualified_email_address_list {

    my ($logged_users) = @_;

    my $users;

    foreach my $user (@{$logged_users}){
	
	if ($user eq 'chado_admin'){
	    next;
	}
	if ($user eq 'access'){
	    next;
	}

	if ($user !~ /\@/){
	    $user .= "\@tigr.org";
	}

	$users.=$user.",";

    }

    # remove the trailing comma
    chop $users;

    return $users;
}


#-----------------------------------------------------------------
# verify_action()
#
#-----------------------------------------------------------------
sub verify_action {

    my ($action) = @_;

    if (!defined($action)){
	$logger->logdie("action was not defined");
    }
    else {
	if (($action eq 'login') || ($action eq 'logout')){
	    if ($logger->is_debug()){
		$logger->debug("action was specified as '$action'");
	    }
	}
	else {
	    $logger->logdie("unexpected action '$action'");
	}
    }

}

#--------------------------------------------------------------
# verify_repository_root()
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



#------------------------------------------------------
# print_usage()
#
#------------------------------------------------------
sub print_usage {

    print STDERR "SAMPLE USAGE:  $0 -D database -U username -P password -a action -c component [-d debug_level] [-h] [-l log4perl] [-m] -r repository -p pipeline_id\n".
    "  -D|--database            = Name of database\n".
    "  -U|--username            = Login name\n".
    "  -P|--password            = Login password\n".
    "  -a|--action              = Action to be taken either 'create' or 'remove'\n".
    "  -c|--component           = Name of workflow component\n".
    "  -d|--debug_level         = Optional - Coati::Logger log4perl logging level.  Default is 0\n".
    "  -h|--help                = Optional - Display pod2usage help screen\n".
    "  -l|--log4perl            = Optional - Log4perl log file (default: /tmp/notify_database_users.pl.log)\n".
    "  -m|--man                 = Optional - Display pod2usage pages for this utility\n".
    "  -w|--repository          = The project's repository root\n".
    "  -p|--pipeline_id         = The pipeline.xml.instance\n";
    exit 1;

}

#------------------------------------------------------
# set_logger()
#
#------------------------------------------------------
sub set_logger {

    my ($debug_level, $log4perl) = @_;

    if (!defined($debug_level)){
	$debug_level = 5;
    }

    if (!defined($log4perl)){
	$log4perl = "/tmp/notify_database_users.pl.log";
    }
    
    my $mylogger = new Workflow::Logger('LOG_FILE'=>$log4perl,
					'LOG_LEVEL'=>$debug_level);
    
    my $logger = Workflow::Logger::get_logger(__PACKAGE__);
    
    return $logger;
}

#------------------------------------------------------
# check_arguments()
#
#------------------------------------------------------
sub check_arguments {

    my ($username, $password, $database, $action, $component, $repository, $pipeline) = @_;

    my $errorctr=0;
    
    if (!defined($username)){
	print STDERR "username was not defined\n";
	$errorctr++;
    }

    if (!defined($password)){
	print STDERR "password was not defined\n";
	$errorctr++;
    }

    if (!defined($database)){
	print STDERR "database was not defined\n";
	$errorctr++;
    }
     
    if (!defined($action)){
	print STDERR "action was not defined\n";
	$errorctr++;
    }
    else {
	#
	# action must be specified
	#
	&verify_action($action);
    }


    if (!defined($component)){
	print STDERR "component was not defined\n";
	$errorctr++;
    }
    
    if (!defined($repository)){
	print STDERR "repository was not defined\n";
	$errorctr++;
    }
    else {
	#
	# repository must be specified
	#
	&verify_repository_root($repository);
    }


    if (!defined($pipeline)){
	print STDERR "pipeline was not defined\n";
	$errorctr++;
    }

   
    if ($errorctr > 0 ) {
	&print_usage();
    }
}
