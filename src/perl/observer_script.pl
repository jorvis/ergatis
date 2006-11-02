#!/usr/local/bin/perl
use lib (@INC,$ENV{"PERL_MOD_DIR"});
no lib "$ENV{PERL_MOD_DIR}/i686-linux";
no lib ".";
#------------------------------------------------------------------------------
#
# author:   Jay Sundaram sundaram@tigr.org
#
# date:     2005-10-27
#
# cvs:      ergatis/workflow/observer_script.pl
#
# $Id$
#
#
#
# Test sample invocation:
#
# perl -I ../lib/ ./observer_script.pl --name=legacy2bsml --event=failure --time=now --id=3434343 --file=/usr/local/scratch/sundaram --props=observer_script.props
#
#
#------------------------------------------------------------------------------
=head1 NAME

observer_script.pl - documentation forthcoming

=head1 SYNOPSIS

USAGE:  observer_script.pl --name <name> --ID <id> --time <time> --event <event> --file <file> --props <props> --message <message> --retval <retval>

=head1 OPTIONS

=over 8

=item B<--name>



=item B<--ID>



=item B<--time>



=item B<--event>



=item B<--props>

    Properties file e.g. ergatis/workflow/observer_script.props

=item B<--message>



=item B<--retval>
 
=back

=head1 DESCRIPTION

    observer_script.pl - documentation forthcoming

    Assumptions:
    1. User has appropriate permissions (to execute script).
    2. All software has been properly installed, all required libraries are accessible.

    Sample usage:
    ./observer_script.pl --name --id --event --time --file --props --message --retval

=head1 CONTACT

    Jay Sundaram
    sundaram@tigr.org

=cut

use Mail::Mailer;
use strict;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use Pod::Usage;
use Config::IniFiles;


$|=1;


#-----------------------------------------------------------------------------
# Eventually, these two values can be read from the properties file
#
my $debug_level = 5;

my $log4perl    = "/usr/local/devel/ANNOTATION/cas/datamanagement/observer/$$.observer_script.pl.log";

#
#---------------------------------------------------------------------------


#-------------------------------------------------------------
# Parse command line options
#-------------------------------------------------------------

my ($name, $id, $time, $event, $file, $props, $message, $retval, $host);


my $results = GetOptions (
			  'name=s'    => \$name,
			  'ID=s'      => \$id, 
			  'time=s'    => \$time,
			  'event=s'   => \$event,
			  'file=s'    => \$file,
			  'props=s'   => \$props,
			  'message=s' => \$message,
			  'host=s'    => \$host,
			  'retval=s'  => \$retval
			  );


print STDERR "name was not specified" if (!defined($name));
print STDERR "event was not specified" if (!defined($event));
print STDERR "time was not specified" if (!defined($time));
print STDERR "file was not specified" if (!defined($file));
print STDERR "props was not specified" if (!defined($props));
print STDERR "ID was not specified" if (!defined($id));



&print_usage() if (!$name or !$event or !$time or !$file or !$props or !$id);

&is_file_readable($props);

my ($config) = &parse_config($props);


my $workflow_component_name = &get_workflow_component_name_from_pipeline_file($file);

my $emaillist = &make_email_list($config, $event, $workflow_component_name);

my $subject = "[$workflow_component_name] [$name] [$event] observer_script.pl";
my $body    = "name:               '$name'\n".
"pipeline:           '$file'\n".
"event:              '$event'\n".
"time:               '$time'\n".
"pipeline:           '$id'\n".
"message:            '$message'\n".
"retval:             '$retval'\n".
"host:               '$host'\n".
"\n\nDocumentation regarding this email notice is available at\n".
"http://intranet/ifx/devel/workflow/user_guide.html\n".
"under the sections titled 'RunWorkflow', 'Observers' and 'Observer Scripts'.\n".
"\n\nContact sundaram\@tigr.org if you have any questions regarding this email notice.";




&send_notification( $emaillist, $subject, $body );




#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#
#                                  END OF MAIN  -- SUBROUTINES FOLLOW
#
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

sub get_workflow_component_name_from_pipeline_file {

    my $file = shift;

    my $workflow_component_name;


    if ($file =~ /\/Workflow\/(\S+)\//){

	$workflow_component_name = $1;
	if ($workflow_component_name =~ /\//){
	    my @tmp = split(/\//, $workflow_component_name);
	    $workflow_component_name = $tmp[0];
	}
    }
    else{
	die "Could not parse the workflow component name from '$file'";
    }

    return ($workflow_component_name);
}


sub parse_config {

    my $props = shift;

    my $cfg = new Config::IniFiles( -file => "$props" );

    my $config = $cfg->{'v'};
    

    return ($config);
}


sub make_email_list {

    my ($config, $event, $name) = @_;


#     testing:   
#    my $outfile = "/usr/local/scratch/sundaram/outfile.txt";
#    open (OUTFILE, ">$outfile") or die "Could not open outfile '$outfile'";
#    print OUTFILE "event '$event' name '$name'\n";


    if (( exists $config->{$name}->{$event}) && (defined($config->{$name}->{$event})) && ($config->{$name}->{$event} ne '0')) {
	return $config->{$name}->{$event};
    }
    else{
	#
	# No need to send notification for this component-event
	#
	exit(0);
    }
}

#------------------------------------------------------
# send_notification()
#
#------------------------------------------------------
sub send_notification {

    my ( $emaillist, $subject, $body ) = @_;

    my $mailer = Mail::Mailer->new ('sendmail');

    $mailer->open({
	             To      => $emaillist,
		     Subject => $subject,
		     From    => undef
		 }) or die("Could not create and send message");
    
    print $mailer $body;
    
    $mailer->close;

}


#-------------------------------------------------------------------
# is_file_readable()
#
#-------------------------------------------------------------------
sub is_file_readable {

    my ( $file) = @_;

    die "file was not defined" if (!defined($file));

    if (!-e $file){
	die "file '$file' does not exist";
    }
    else{
	if ((!-r $file)){
	    die "file '$file' does not have read permissions";
	}
	if ((-z $file)){
	    die "file '$file' has no content";
	}
    }

    return 1;
   

}#end sub is_file_readable()

#------------------------------------------------------
# print_usage()
#
#------------------------------------------------------
sub print_usage {

    print STDERR "SAMPLE USAGE:  $0 --name <name> --ID <id> --time <time> --event <event> --file <file> --props <props> --message <message> --retval <retval>\n".
    "  --name         = name\n".
    "  --ID           = id\n".
    "  --time         = time\n".
    "  --event        = event\n".
    "  --file         = file\n".
    "  --props        = props\n".
    "  --message      = message\n".
    "  --host         = host\n".
    "  --retval       = retval\n";
    exit 1;

}

