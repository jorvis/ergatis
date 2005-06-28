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
# usage:   CAS group's run_wf.sh script should configure and deploy this script as part of the
#          RunWorkflow invocation (--scripts command-line argument)
#
#
#
#-----------------------------------------------------------------------------------------------

use strict;
use Mail::Mailer;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use Pod::Usage;
use Config::IniFiles;
use Data::Dumper;
$|=1;

umask(0000);

#-------------------------------------------------------------
# Parse command line options
#-------------------------------------------------------------

my ($name, $id, $time, $event, $file, $props, $sample);

my $results = GetOptions (
			  'name=s'        => \$name,
			  'id=s'          => \$id,
			  'time=s'        => \$time,
			  'event=s'       => \$event,
			  'file=s'        => \$file,
			  'props=s'       => \$props,
			  'sample'        => \$sample,
			  );

print STDERR "name was not defined\n"  if (!defined($name));
print STDERR "id was not defined\n"    if (!defined($id));
print STDERR "time was not defined\n"  if (!defined($time));
print STDERR "event was not defined\n" if (!defined($event));
print STDERR "file was not defined\n"  if (!defined($file));
print STDERR "props was not defined\n" if (!defined($props));

&print_usage() if (!$name or !$id or !$time or !$event or !$file or !$props);
&print_sample() if $sample;


#die "name '$name' id '$id' time '$time' event '$event' file '$file' props '$props'";

my $addiprops = &get_additional_properties($props);
#print Dumper $addiprops;die;



my $nameevent = $name . '_' . $event;
#die "nameevent:$nameevent";

#
# default subject and body
#
my $subject = "[Workflow Observer] id '$id' name '$name' ";
my $body = "name '$name'\nid '$id'\ntime '$time'\nevent '$event'\nfile '$file'\n";


#
# append additional event specific information to the subject and body
#
if ((exists $addiprops->{'v'}->{$nameevent}->{'subject'}) and (defined($addiprops->{'v'}->{$nameevent}->{'subject'}))){
    $subject .= $addiprops->{'v'}->{$nameevent}->{'subject'};
}
if ((exists $addiprops->{'v'}->{$nameevent}->{'body'}) and (defined($addiprops->{'v'}->{$nameevent}->{'body'}))){
    $body .= $addiprops->{'v'}->{$nameevent}->{'body'};
}

my $final_to = {};

#
# Collect all the default to email addresses in a hash so that we can ensure uniqueness
#
if ((exists $addiprops->{'v'}->{'to_list'}->{'to'}) and (defined($addiprops->{'v'}->{'to_list'}->{'to'}))){
    
    my $to_list = $addiprops->{'v'}->{'to_list'}->{'to'};
    my @sendlist = split(/,/,$to_list);
    foreach my $name (sort @sendlist){
	$final_to->{$name} = 1;
    }
}

#
# Now collect all name-event specific address in the same hash
#
if ((exists $addiprops->{'v'}->{$nameevent}->{'to'}) and (defined($addiprops->{'v'}->{$nameevent}->{'to'}))){

    my $to_list = $addiprops->{'v'}->{$nameevent}->{'to'};
    my @sendlist = split(/,/,$to_list);
    foreach my $name (sort @sendlist){
	if (!exists $final_to->{$name}){
	    $final_to->{$name} = 1;
	}
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

#print STDERR ("subject: $subject\n");
#print STDERR ("body:    $body\n");
#print STDERR ("to_list: $to_list\n");
#die;


&send_notification($to_list, $subject, $body);
exit(0);




#----------------------------------------------------------------------
# send_notification()
#
#----------------------------------------------------------------------
sub send_notification {

    my ($to_list, $subject, $body) = @_;

    my $mailer = Mail::Mailer->new ('mail');
    $mailer->open({
	To   => "$to_list",
	Subject => "Workflow $$ completed"
    }) or die "Could not create and send message";
    
    print $mailer $body;
    
    $mailer->close;

    print "Notification sent\n";

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


#-----------------------------------------------------------------
# print_sample()
#
#-----------------------------------------------------------------
sub print_sample {

    print<<END_SAMPLE;

[to_list]
to=sundaram\@tigr.org

[legacy2chado_start]
body=
subject=
to=sundaram\@tigr.org,jay_sundaram\@hotmail.com

[legacy2chado_finish]
body=
subject=
to=

[legacy2chado_failure]
body=
subject=
to=jay_sundaram\@hotmail.com,2409947066\@vtext.com,sundaram\@tigr.org

[legacy2chado_interrupt]
body=
subject=

[legacy2chado_resume]
body=
subject=

[legacy2chado_restart]
body=
subject=


[chado2bsml_start]
body=
subject=

[chado2bsml_finish]
body=
subject=

[chado2bsml_failure]
body=
subject=

[chado2bsml_interrupt]
body=
subject=

[chado2bsml_resume]
body=
subject=

[chado2bsml_restart]
body=
subject=

END_SAMPLE


exit(0);

}


#-----------------------------------------------------------------
# print_usage()
#
#-----------------------------------------------------------------
sub print_usage {

    print STDERR "SAMPLE USAGE:  $0 --name=<name> --id=<id> --time=<time> --event=<event> --file=<file> --props=<props> [--sample]\n";
    print STDERR "  --name     = \n";
    print STDERR "  --id       = \n";
    print STDERR "  --time     = \n";
    print STDERR "  --event    = \n";
    print STDERR "  --file     = \n";
    print STDERR "  --props    = \n";
    print STDERR "  --sample   = Optional - to view sample props file\n";
    exit 1;
   


}
