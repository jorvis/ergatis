#!/usr/local/bin/perl
#print STDERR "Starting top\n";
use strict;

use CGI qw(:standard);

use Tree::DAG_Node;
use File::Find;
use File::Basename;

my $instancexml = param('instancexml');
my $inifile = param('inifile');
my $template = param('xmltemplate');

my $log =  $instancexml.".log";
my $out = $instancexml.".out";

my $child_pid;
$ENV{'HTTP_USER_AGENT'} = '';
if(($inifile ne "") && ($template ne "")){
    my $child_pid = fork;
    if($child_pid){
	while(! (-e $instancexml)){
	    sleep 3;
	}
	print redirect(-uri=>"./show_pipeline.cgi?xmltemplate=$instancexml");    
	exit;
    }
    else{
	close STDOUT;
	close STDERR;
	close STDIN;
        #  Fork again.  This helps separate the background process from
        #  the httpd process.  If we're in the original child, $gpid will
        #  hold the process id of the "grandchild", and if we're in the
        #  grandchild it will be zero.
	my $gpid = fork;
	if(! $gpid){
            #  We're in the grandchild.
	    chdir '/usr/local/scratch/workflow';
	    use POSIX qw(setsid);
	    setsid() or die "Can't start a new session: $!";
	    $ENV{'WF_ROOT'} = "/usr/local/devel/ANNOTATION/workflow-2.2B2";
	    $ENV{'WF_ROOT_INSTALL'} = "/usr/local/devel/ANNOTATION/workflow-2.2B2";
	    $ENV{'LD_LIBRARY_PATH'} = "/usr/local/lib:/usr/local/packages/sybase/OCS/lib:/usr/local/packages/sybase/lib";
	    $ENV{'SYBASE'} = "/usr/local/packages/sybase";
	    $ENV{'PATH'} = "$ENV{'WF_ROOT'}:$ENV{'WF_ROOT'}/bin:$ENV{'WF_ROOT'}/add-ons/bin:$ENV{'PATH'}";
	    my $createexecstr = "$ENV{'WF_ROOT'}/CreateWorkflow -debug -i $instancexml -t $template -c $inifile --delayedbuild=true --autobuild=false > $instancexml.create.out 2>&1";
	    system("$createexecstr");
	    my $runexecstr = "$ENV{'WF_ROOT'}/RunWorkflow -i $instancexml > $instancexml.run.out 2>&1";
	    system($runexecstr);
	    exit;
	}
	exit;
    }
}
else{
    my $child_pid = fork;
    if($child_pid){
	print redirect(-uri=>"./show_pipeline.cgi?xmltemplate=$instancexml");    
	exit;
    }
    else{
	close STDOUT;
	close STDERR;
	close STDIN;
        #  Fork again.  This helps separate the background process from
        #  the httpd process.  If we're in the original child, $gpid will
        #  hold the process id of the "grandchild", and if we're in the
        #  grandchild it will be zero.
	my $gpid = fork;
	if(! $gpid){
	    chdir '/usr/local/scratch/workflow';
	    use POSIX qw(setsid);
	    setsid() or die "Can't start a new session: $!";

	    $ENV{'WF_ROOT'} = "/usr/local/devel/ANNOTATION/workflow-2.2B2";
	    $ENV{'WF_ROOT_INSTALL'} = "/usr/local/devel/ANNOTATION/workflow-2.2B2";
	    $ENV{'SYBASE'} = "/usr/local/packages/sybase";
	    $ENV{'PATH'} = "$ENV{'WF_ROOT'}:$ENV{'WF_ROOT'}/bin:$ENV{'WF_ROOT'}/add-ons/bin:$ENV{'PATH'}";
	    $ENV{'LD_LIBRARY_PATH'} = "/usr/local/lib:/usr/local/packages/sybase/OCS/lib:/usr/local/packages/sybase/lib";
	    my $runstring = "$ENV{'WF_ROOT'}/RunWorkflow -i $instancexml -l $instancexml.log -v5 >& $instancexml.run.out";
	    my $rc = 0xffff & system($runstring);
	    use Data::Dumper;
	    print STDERR Dumper(%ENV);
	    printf STDERR "system(%s) returned %#04x: $rc ","$runstring", $rc;
	    if($rc == 0){
		print STDERR "ran with normal exit\n";
	    }
	    elsif( $rc == 0xff00){
		print STDERR "command failed: $!\n";
	    }
	    elsif(($rc & 0xff) == 0) {
		$rc >>= 8;
		print STDERR "ran with no-zero exit status $rc\n";
	    }
	    else {
		print STDERR "ran with ";
		if($rc & 0x80){
		    $rc &= ~0x80;
		    print STDERR "coredump from ";
		}
		print STDERR "signal $rc\n";
	    }
	}
	exit;
    }
}

exit;

#use lib "/home/condor/production/request/lib";
#use TIGR::HTCRequest;


#print "hello\n";

#my $request = TIGR::HTCRequest->new(group => "antware",
#				    initialdir => "/usr/local/annotation/FIXGHI/BER2multi_dir/WORKING");
#$request->username('angiuoli');
#$request->set_command("/home/dsommer/working/test.sh");
#$dd = $request->get_command();
#print " command is $dd \n";
#$request->add_param("-output=GCG");
#$request->add_param("-infile", "/home/dsommer/tmp/.*.fa\$", DIR);
#$request->add_param("-outfile", "/usr/local/annotation/FIXGHI/BER2multi_dir/WORKING/\$(Name).msf");

#$request->set_output("/home/dsommer/working/htc_test.\$(Name).out");
#$request->set_error("/home/dsommer/working/htc_test.\$(Name).err");


#$request->add_param("filekey", "/home/dsommer/working/filelines", FILE);


#$request->set_getenv(1);
#$request->set_env_list(...);

#my $xx = $request->to_xml();
#print $xx;

#my $id = $request->submit();
#print " Request id was $id \n";

# wait for job finish
#$request->wait_for_request();

#my $message = $request->get_message();
#if($request->get_state() eq "FAILURE") {
#    print " dpsSearch failed: $message \n";
#    exit(1);
#} elsif ($request->get_state() eq "FINISHED") {
#    print " Request finished correctly\n";
#} else {
#    print " Request finished with state $request->get_state() and $message \n";
#}

#print "command finished\n";

