#!/usr/bin/perl
BEGIN{foreach (@INC) {s/\/usr\/local\/packages/\/local\/platform/}};
use lib (@INC,$ENV{"PERL_MOD_DIR"});
no lib "$ENV{PERL_MOD_DIR}/i686-linux";
no lib ".";

use strict;
BEGIN {
use Workflow::Repository;
}
use Date::Manip;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);

my %options;
my $results = GetOptions (\%options, 
                          'directory|d=s', 
                          'type|t=s', 
			  'log|l=s',
                          'debug=s', 
			  'skiprun',
                          'help|h' );

my $logfile = $options{'log'} || Ergatis::Logger::get_default_logfilename();

my $logger = new Ergatis::Logger('LOG_FILE'=>$logfile,
				  'LOG_LEVEL'=>$options{'debug'});

# display documentation
if( $options{'help'} ){
    pod2usage( {-exitval=>0, -verbose => 2, -output => \*STDOUT} );
}

my $wfrepository  = new Workflow::Repository('REPOSITORY_ROOT'=> $options{'directory'});

my @workflowtypes = $wfrepository->get_types_dirs();
    
foreach my $type (@workflowtypes){
    if($type eq $options{'type'} || (!$options{'type'})){
	my @workflowinstances = $wfrepository->get_instance_dirs($type);
	foreach my $instance (@workflowinstances){
	    my $instancexml = $wfrepository->get_instance_xml($type,$instance);
	    
	    if(-e $instancexml){
		print "Type: $type \n";
		print "\tUID: $instance $instancexml\n";
		
		my($topstate,$topstarttime,$topendtime,$stateslookup,$timeslookup,$hostslookup) = $wfrepository->get_workflow_status($instancexml);
		print "\tState: $topstate\n";
		print "\tStartTime: ",&UnixDate(&DateCalc("Jan 1, 1970  00:00:00 GMT",$topstarttime),"%c"),"\n";
		print "\tEndTime: ",&UnixDate(&DateCalc("Jan 1, 1970  00:00:00 GMT",$topendtime),"%c"),"\n";

		my $currtime = ($topendtime ne "") ? $topendtime : &UnixDate("today","%s");
		my $elapsedtime = $currtime - $topstarttime;
						
		print "\tElapsedTime: ",&Delta_Format(&ParseDateDelta($elapsedtime),1,"%ht"),"hrs\n";

		my $totalstates;
		foreach my $state (keys %$stateslookup){
		    $totalstates += $stateslookup->{$state};
		}
		foreach my $state (keys %$stateslookup){
		    print "\t\tState $state: $stateslookup->{$state}/$totalstates (";
		    printf("%.1f", ($stateslookup->{$state}/$totalstates)*100);
		    print "%)\n";
		}
		foreach my $times (keys %$timeslookup){
		    my $timeperjob = sprintf("%d",($timeslookup->{$times}/$stateslookup->{'complete'}));
		    print "\t\tTime $times: ",&Delta_Format(&ParseDateDelta("$timeslookup->{$times}"),1,"%ht"),"hrs (",&Delta_Format(&ParseDateDelta("$timeperjob"),1,"%mt"),"min/job)\n";
		}
		
		my $speedup = ($timeslookup->{'exectime'}/$elapsedtime);
		
		print "\tSpeedup (exectime/elapsedtime) = ";
		printf("%.1f", $timeslookup->{'exectime'}/$elapsedtime);
		print "x\n";

		print "\tNumber of hosts used ",scalar(keys %$hostslookup)+1,"\n";
	    }
	}
    }
}


	    
 
