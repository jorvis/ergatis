package Workflow::Repository;

# $Id$

# Copyright (c) 2002, The Institute for Genomic Research. All rights reserved.

=head1 NAME

Repository.pm - A module for managing a Workflow repository

=head1 VERSION

This document refers to version $Name$ of frontend.cgi, $Revision$. 
Last modified on $Date$

=head1 SYNOPSIS

=head1 DESCRIPTION

my $workflowrepository = new Workflow::Repository('PATH'=>$repositorypath);

=over 4

=cut


use strict;
use Data::Dumper;
BEGIN {
    require '/usr/local/devel/ANNOTATION/ard/chado-v1r5b1/lib/site_perl/5.8.5/Workflow/Logger.pm';
    import Workflow::Logger;
}
use Date::Manip;
use XML::Twig;
use Data::Dumper;

=item new

B<Description:> The module constructor.

B<Parameters:> %arg, a hash containing attribute-value pairs to
initialize the object with. Initialization actually occurs in the
private _init method.

my $bsmlrepository = new Workflow::Repository('PATH'=>$repositorypath);

B<Returns:> $self (A Workflow::Repository object).

=cut

sub new {
    my ($class) = shift;
    my $self = bless {}, ref($class) || $class;
    $self->{_logger} = Workflow::Logger::get_logger(__PACKAGE__);
    $self->{_WORKFLOW_SUBDIR} = "Workflow";
    $self->{_INSTANCE_FILENAME} = "master.xml";
    $self->_init(@_);
    $self->{"_PATH"} = $self->{"_REPOSITORY_ROOT"}."/".$self->{_WORKFLOW_SUBDIR}."/";
    $self->{_logger}->debug("Setting repository path $self->{_PATH}") if($self->{_logger}->is_debug());
    $self->_create_directory($self->{"_PATH"});
    return $self;
}


=item $obj->_init([%arg])

B<Description:> Initialize object variables

B<Parameters:> %arg, a hash containing attributes to initialize the testing
object with. Keys in %arg will create object attributes with the same name,
but with a prepended underscore.

B<Returns:> None.

=cut

sub _init {
    my $self = shift;
    
    my %arg = @_;
    $self->{_logger}->debug(Dumper(@_)) if($self->{_logger}->is_debug());
    foreach my $key (keys %arg) {
	$self->{_logger}->debug("Parsing argument $key=$arg{$key}") if($self->{_logger}->is_debug());
        $self->{"_$key"} = $arg{$key};
    }
    if(!($self->{"_REPOSITORY_ROOT"})){
	$self->{_logger}->logdie("Required parameter REPOSITORY_ROOT not passed to object constructor");
    }
}

#return the name of the Workflow repository
sub get_root_dir{
    my $self = shift;
    return $self->{"_PATH"};
}

#return all type dirs
sub get_types_dirs{
    my $self = shift;
    opendir WORKFLOWDIRS, "$self->{_PATH}" or $self->{_logger}->logdie("Can't read directory $self->{_PATH}");
    my @workflowtypes =  grep /\w+$/,readdir WORKFLOWDIRS; 
    closedir WORKFLOWDIRS;
    return @workflowtypes;
}

sub get_instance_dirs{
    my $self = shift;
    my $type = shift;
    opendir WORKFLOWSUBDIRS, "$self->{_PATH}/$type" or $self->{_logger}->logdie("Can't read directory $self->{_PATH}/$type");
    my @workflowinstances = grep /\d+$/,readdir WORKFLOWSUBDIRS;
    closedir WORKFLOWDIRS;
    return @workflowinstances;
}

sub get_instance_xml{
    my $self = shift;
    my $type = shift;
    my $instanceid = shift;
    
    return "$self->{_PATH}/$type/$instanceid/$self->{_INSTANCE_FILENAME}";
}

sub create_working_dir{
    my $self = shift;
    my $type = shift;
    my $uid = shift;
    $type =~ s/\s//g;
    $uid =~ s/\s//g;
    my $dir = "$self->{_PATH}/$type/$uid";
    $self->_create_directory($dir);
    return $dir;
}

sub get_all_instances{
    my $self = shift;
    my @workflows;

    my @workflowtypes = $self->get_types_dirs();
    
    foreach my $type (@workflowtypes){
	my @workflowinstances = $self->get_instance_dirs($type);
	foreach my $instance (@workflowinstances){
	    return $self->_parse_workflow_instance_dir($self->{_PATH}/$type/$instance);
	    opendir WORKFLOWINSTANCEDIRS, "$self->{_PATH}/$type/$instance" or $self->{_logger}->logdie("Can't read directory $self->{_PATH}/$type/$instance");
	    my $instancefile = "$self->{_PATH}/$type/$self->{_INSTANCE_FILENAME}";
	    my $directory = "$self->{_PATH}/$type";
	    if(! (-e $instancefile)){
		 $self->{_logger}->logdie("Can't read instance file $instancefile");
	     }
	    my $status = &get_workflow_status($instancefile);
	    my $states = &get_workflow_states($instancefile);
	    my ($starttime,$endtime) = &get_workflow_runtimes($instancefile);

	    push @workflows,$self->_parse_workflow_instance_dir($instancefile);
	}
    }
    return \@workflows;
}

sub get_instance_by_type{
    my $self = shift;
    my $type = shift;
    my @workflows;
    opendir WORKFLOWDIRS, "$self->{_PATH}" or $self->{_logger}->logdie("Can't read directory $self->{_PATH}");
    my @workflowtypes = readdir WORKFLOWDIRS; 
    foreach my $type (@workflowtypes){
	opendir WORKFLOWSUBDIRS, "$self->{_PATH}/$type" or $self->{_logger}->logdie("Can't read directory $self->{_PATH}/$type");
	my @workflowinstances = readdir WORKFLOWSUBDIRS;
	foreach my $instancefile (@workflowinstances){
	     push @workflows,$self->_parse_workflow_instance_dir($instancefile);
	}
    }
    return \@workflows;
}

sub get_instance_by_command_id{
#search for workflow instance containing command
#can we return all parent instances as well?

}

sub get_component_list{
    my ($self,$dir) = @_;
    opendir COMPONENTDIRS, "$dir" or $self->{_logger}->logdie("Can't read directory $dir");
    my @components =  grep /conf$/,readdir COMPONENTDIRS;
    return @components;
}

sub _parse_workflow_instance_dir{
    my $self = shift;
    my $instance_dir = shift;
    my $type = shift;
    
    opendir WORKFLOWINSTANCEDIRS, "$instance_dir" or $self->{_logger}->logdie("Can't read directory $instance_dir");
    my $instancefile = "$instance_dir/$self->{_INSTANCE_FILENAME}";
    if(! (-e $instancefile)){
	$self->{_logger}->logdie("Can't read instance file $instancefile");
    }
    my $status = &get_workflow_status($instancefile);
    my $states = &get_workflow_states($instancefile);
    my ($starttime,$endtime) = &get_workflow_runtimes($instancefile);
    
    my $workflowref = {'name'=>$type,
		       'filename'=>$instancefile,
		       'directory'=>$instance_dir,
		       'state'=>$status,
		       'start_time'=>$starttime,
		       'end_time'=>$endtime,
		       'job_states'=>$states};
    return $workflowref;
}

sub get_workflow_status{
    my($self,$instancefile) = @_;
    $self->{_logger}->debug("Getting status of $instancefile") if($self->{_logger}->is_debug());
    my $topstate;
    my $topstarttime;
    my $topendtime;

    my $t1 = new XML::Twig( TwigHandlers => { 'commandSetRoot/commandSet' =>
						  sub {
						      my ($t, $elt) = @_;
						      
						      $topstate = $elt->first_child('state')->text();
						      $self->{_logger}->debug("Top state $topstate") if($self->{_logger}->is_debug());
						      
						      $topstarttime = &UnixDate($elt->first_child('startTime')->text(),"%s");
						      $self->{_logger}->debug("Starttime $topstarttime") if($self->{_logger}->is_debug());
						      
						      if($elt->first_child('endTime')){
							  $topendtime = &UnixDate($elt->first_child('endTime')->text(),"%s");
							  $self->{_logger}->debug("Endtime $topendtime") if($self->{_logger}->is_debug());
						      }
						      return 1;
						  }}
			    );


    my @filenames = "$instancefile";

    $self->{_logger}->debug("Parsing $instancefile") if($self->{_logger}->is_debug());

    $t1->parsefile("$instancefile");


    my $stateslookup = {};
    my $executionhosts = {};
    my $timeslookup = {};

    my $uniquefiles = {};

    foreach my $filename (@filenames){
	if(-e $filename && (! exists($uniquefiles->{$filename}))){
	    $self->{_logger}->debug("Parsing referenced xml file $filename for filename references") if($self->{_logger}->is_debug());
	    
	    #add nested filename refs to stack
	    my $t3 = new XML::Twig( TwigRoots => {'fileName' => 
							 sub {
							     my ($t, $elt) = @_;
							     push @filenames, $elt->text();
							     return 1;
							 }}
				    );
	    
	    eval{
		$t3->parsefile($filename);
	    };

	    $uniquefiles->{$filename} = 1;
	}
    }

    foreach my $filename (keys %$uniquefiles){
	#queueing time
	my $totalqueuetime = 0;
	#est. exectime
	my $totalexectime = 0;
	#est. distributed overhead
	my $totaldoverheadtime = 0;
	
	my $t2= new XML::Twig( TwigRoots => { 'command' => 
						 sub {
						     my ($t, $elt) = @_;
						     
						     
						     my ($state) = $elt->first_child('state')->text();
						     $stateslookup->{$state}++;
						     
						
						     my ($starttime,$endtime) = _get_command_times($elt);
						     if($state eq "complete" && ($endtime != 0)){
							 $self->{_logger}->debug("Command state $state starttime:$starttime endtime:$endtime") if($self->{_logger}->is_debug());
							 my ($dce) = $elt->first_child('dceSpec');
							 if($dce){
							     my $host = $dce->first_child('executionHost')->text();
							     $self->{_logger}->debug("executionHost $host") if($self->{_logger}->is_debug());
							     $executionhosts->{$host}++;
							     my $log = $dce->first_child('log')->text();
							     $self->{_logger}->debug("condor log $log") if($self->{_logger}->is_debug());
							     my ($submittime,$executionstart,$executionend) = _parse_condor_logs($log);
							     $self->{_logger}->debug("submittime:$submittime executionstart:$executionstart executionend:$executionend") if($self->{_logger}->is_debug());
							     my $exectime = ($executionend-$executionstart);
							     if($exectime < 0){
								 print STDERR "Possible clock sync problem executionstart ($executionstart):", 
								 &UnixDate(&DateCalc("Jan 1, 1970  00:00:00 GMT",$executionstart)), 
								 " executionend ($executionend):",
								 &UnixDate(&DateCalc("Jan 1, 1970  00:00:00 GMT",$executionend)),
								 " ",$elt->first_child('id')->text(),"\n";
							     }
							     else{
								 $totalexectime += $exectime;
							     }
							     
							     my $doverheadtime1 = ($submittime - $starttime);
							     if($doverheadtime1 < 0){
								 print STDERR "Possible clock sync problem startTime ($starttime):", 
								 &UnixDate(&DateCalc("Jan 1, 1970  00:00:00 GMT",$starttime)), 
								 " submit time ($submittime):",
								 &UnixDate(&DateCalc("Jan 1, 1970  00:00:00 GMT",$submittime)),
								 " $elt->first_child('id')->text();\n";
							     }
							     else{
								 my $doverheadtime2 = ($endtime-$executionend);
								 if($doverheadtime2 < 0){
								     print STDERR "Possible clock sync problem endTime ($endtime):", 
								     &UnixDate(&DateCalc("Jan 1, 1970  00:00:00 GMT",$endtime)), 
								     " execution end time ($executionend):",
								     &UnixDate(&DateCalc("Jan 1, 1970  00:00:00 GMT",$executionend)),
								     " ",$elt->first_child('id')->text(),"\n";
								 }else{
								     $totaldoverheadtime += ($doverheadtime1 + $doverheadtime2);
								 }
							     }
							     
							     my $queuetime = ($executionstart - $submittime);
							     if($queuetime < 0){
								 print STDERR "Possible clock sync problem submittime ($submittime)", 
								 &UnixDate(&DateCalc("Jan 1, 1970  00:00:00 GMT",$submittime)), 
								 " execution start time ($executionstart):",
								 &UnixDate(&DateCalc("Jan 1, 1970  00:00:00 GMT",$executionstart)),
								 " ",$elt->first_child('id')->text(),"\n";
							     }else{
								 $totalqueuetime += $queuetime;
							     }
							     $self->{_logger}->debug("totalexectime:$totalexectime executionstart:$totalexectime executionend:$totalqueuetime") if($self->{_logger}->is_debug());
							 }
							 else{
							     my $exectime = ($endtime-$starttime);
							     if($exectime < 0){
								 print STDERR "Possible clock sync problem starttime ($starttime):", 
								 &UnixDate(&DateCalc("Jan 1, 1970  00:00:00 GMT",$starttime)), 
								 " endtime ($endtime):",
								 &UnixDate(&DateCalc("Jan 1, 1970  00:00:00 GMT",$endtime)),
								 " ",$elt->first_child('id')->text(),"\n";
							     }else{
								 $totalexectime += ($endtime-$starttime);
							     }
							 }
						     }
						     return 1;
						 }
					  }
			       );
			       
	eval{
	    $t2->parsefile($filename);
	};
	
	$timeslookup->{'exectime'} += $totalexectime;
	$timeslookup->{'queuetime'} += $totalqueuetime;
	$timeslookup->{'distributed_overhead'} += $totaldoverheadtime;
	
    }
    return ($topstate,$topstarttime,$topendtime,$stateslookup,$timeslookup,$executionhosts);
}

sub _get_command_times{
    my( $elt)= @_;
    my $starttime= 0;
    my $endtime=0;
    if($elt->first_child('startTime')){
	$starttime = $elt->first_child('startTime')->text();
    }
    if($elt->first_child('endTime')){
	$endtime = $elt->first_child('endTime')->text();
    }
    return (&UnixDate($starttime,"%s"),&UnixDate($endtime,"%s"));
}

sub _parse_condor_logs{
    my($log) = @_;
    my $submittime = 0;
    my $executionstart = 0;
    my $executionend = 0;

    open LOG, "$log" or warn("Can't open log file $log for parsing $?");
    while (my $line = <LOG>){
	if($line =~ /Job submitted/){
	    ($submittime) = ($line =~ /([\d\/]+\s+[\d\:]+)\s+Job submitted/);
	}
	elsif($line =~ /Job executing/){
	    ($executionstart) = ($line =~ /([\d\/]+\s+[\d\:]+)\s+Job executing/);
	}
	elsif($line =~ /Job terminated/){
	    ($executionend) = ($line =~ /([\d\/]+\s+[\d\:]+)\s+Job terminated/);
	}
    }
    close LOG;
    return (&UnixDate($submittime,"%s"),&UnixDate($executionstart,"%s"),&UnixDate($executionend,"%s"));
}
    
sub _create_directory{
    my($self, $dir) = @_;
    my $ret = system("mkdir -m 777 -p $dir");
    $ret >>= 8;
    return $ret;
}

return 1;
