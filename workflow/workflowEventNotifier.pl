#!/usr/local/bin/perl -w
#-----------------------------------------------------------------------------------------------
#
# script: workflowEventNotifier.pl 
# author: sundaram
# date:   2004-05-14
#
# purpose: Configurable observer script which receives Workflow event notifications and
#          sends email notification to appropriate individuals to report workflow status
#
#          See: Antware's "Observer Scripts" http://intranet/ifx/devel/workflow/main_frame.html
#
#
# usage:   CAS group's run_wf.sh script should configure and deploy this script
#        
#-----------------------------------------------------------------------------------------------

use strict;
use Mail::Mailer;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use Pod::Usage;
use Config::IniFiles;
use Data::Dumper;
$|=1;

#-------------------------------------------------------------
# Parse command line options
#-------------------------------------------------------------

my ($name, $id, $time, $event, $file, $props);

my $results = GetOptions (
			  'name=s'        => \$name,
			  'id=s'          => \$id,
			  'time=s'        => \$time,
			  'event=s'       => \$event,
			  'file=s'        => \$file,
			  'props=s'       => \$props,
			  );

print STDERR "name was not defined\n"  if (!defined($name));
print STDERR "id was not defined\n"    if (!defined($id));
print STDERR "time was not defined\n"  if (!defined($time));
print STDERR "event was not defined\n" if (!defined($event));
print STDERR "file was not defined\n"  if (!defined($file));
print STDERR "props was not defined\n" if (!defined($props));

&print_usage() if (!$name or !$id or !$time or !$event or !$file or !$props);



#die "name '$name' id '$id' time '$time' event '$event' file '$file' props '$props'";

my $addiprops = &get_additional_properties($props);
print Dumper $addiprops;die;


#
# Notification will only be sent if event is a qualified event
#
foreach my $q_event (sort @{$addiprops->{$name}->{'qualified_events'}}){

    if ($event eq $q_event){
	
	#
	# default subject and body
	#
	my $subject = "[Workflow Observer] id '$id' name '$name' ";
	my $body = "name '$name'\nid '$id'\ntime '$time'\nevent '$event'\nfile '$file'\n";


	#
	# append additional event specific information to the subject and body
	#
	if ((exists $addiprops->{$name}->{$event}->{'subject'}) and (defined($addiprops->{$name}->{$event}->{'subject'}))){
	    $subject .= $addiprops->{$name}->{$event}->{'subject'};
	}
	if ((exists $addiprops->{$name}->{$event}->{'body'}) and (defined($addiprops->{$name}->{$event}->{'body'}))){
	    $body .= $addiprops->{$name}->{$event}->{'body'};
	}



	my $final_to = {};
	
	#
	# check if there are event specific individuals to notify
	#
	if ((exists $addiprops->{$name}->{$event}->{'to'}) and (defined($addiprops->{$name}->{$event}->{'to'}))){

	    foreach my $spec (sort @{$addiprops->{$name}->{$event}->{'to'}}){
		$final_to->{$spec} = '';
	    }
	}

	#
	# build unique 'to list'
	#
	foreach my $to (sort @{$addiprops->{'to_list'}}){
	    if (!exists $final_to->{$to}){
		$final_to->{$to} = '';
	    }
	}

	#
	# generate the to_list
	#
	my $to_list;
	foreach my $key (sort keys %$final_to){
	    $to_list .= $key . ',';
	}
	chop $to_list;

	print STDERR ("to_list: $to_list\n");
	print STDERR ("subject: $subject\n");
	print STDERR ("body:    $body\n");


	&send_notification($to_list, $subject, $body);
	exit(0);
    }
}


#----------------------------------------------------------------------
# send_notification()
#
#----------------------------------------------------------------------
sub send_notification {

    my ($to_list, $subject, $body) = @_;

    my $mailer = Mail::Mailer->new ('mail');
    $mailer->open({
	To   => '$to_list',
	Subject => "Workflow $$ completed"
    }) or die "Could not create and send message";
    
    print $mailer $body;
    
    $mailer->close;
}


#----------------------------------------------------------------------
# get_additional_properties()
#
#----------------------------------------------------------------------
sub get_additional_properties {

    my $file = shift;

    die "'get_additional_properties' file was not defined" if (!defined($file));

    my $hash = new Config::IniFiles( -file => $file );

    die "'get_additional_properties' hash was not defined" if (!defined($hash));

    return $hash;
 
}#end sub get_additional_properties()





# $addiprops = {
#                qualified_events => [
#                                      'start',
#                                      'finish',
#                                      'failure'
#                                     ],
#                to_list => [
#                             'sundaram@tigr.org',
#                             '2409947066@vtext.com',
#                             'angiuoli@tigr.org'
#                           ],
#
#                start     => {
#                                body => 'some default message when start event occurs',
#                                subject => 'some default subject when start event occurs',
#                                to => [
#                                        'angiuoli@tigr.org'
#                                      ]
#                             },
#                finish    => {
#                                body => 'some default message when finish event occurs',
#                                subject => 'some default subject when finish event occurs',
#                                to => [
#                                         'crabtree@tigr.org'
#                                      ]
#                             },
#                failure   => {
#                                body => 'some default message when failure event occurs',
#                                subject => 'some default subject when failure event occurs',
#                                to => [
#                                         'sundaram@tigr.org,
#                                          angiuoli@tigr.org'
#                                      ]
#                             },
#                resume    => {
#                                body => 'some default message when resume event occurs',
#                                subject => 'some default subject when resume event occurs'
#                             },
#                interrupt => {
#                                body => 'some default message when interrupt event occurs',
#                                subject => 'some default subject when interrupt event occurs'
#                             },
#                restart => {
#                                body => 'some default message when restart event occurs',
#                                subject => 'some default subject when restart event occurs'
#                             },
#                subject   => {
#                                body => 'some default message when subject event occurs',
#                                subject => 'some default subject when subject event occurs'
#                             },
#                body      => {
#                                body => 'some default message when body event occurs',
#                                subject => 'some default subject when body event occurs'
#                             }
#              }
#
#
#
# Sample ini props file:
#
#
# [qualified_events]
# start=1
# finish=1
# failure=1
#
# [to_list]
# sundaram@tigr.org=1
# angiuoli@tigr.org=1
# 2409947066@vtext.com=1
#
# [start]
# body=
# subject=
#
# [finish]
# body=
# subject=
#
# [failure]
# body=
# subject=
#
# [interrupt]
# body=
# subject=
#
# [resume]
# body=
# subject=
#
# [restart]
# body=
# subject=
#
#


#-----------------------------------------------------------------
# print_usage()
#
#-----------------------------------------------------------------
sub print_usage {

    print STDERR "SAMPLE USAGE:  $0 --name=<name> --id=<id> --time=<time> --event=<event> --file=<file> --props=<props>\n";
    print STDERR "  --name     = \n";
    print STDERR "  --id       = \n";
    print STDERR "  --time     = \n";
    print STDERR "  --event    = \n";
    print STDERR "  --file     = \n";
    print STDERR "  --props    = \n";
    exit 1;
   


}
