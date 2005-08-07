#!/usr/local/bin/perl
use strict;

use CGI qw(:standard);

use Tree::DAG_Node;
use File::Find;
use File::Basename;

## toggles debugging messages
my $debugging = 0;
my $debug_file = '/tmp/cram.debug';

my $instancexml = param('instancexml');
my $inifile = param('inifile');
my $template = param('xmltemplate');

my $log =  $instancexml.".log";
my $out = $instancexml.".out";

my $child_pid;
$ENV{'HTTP_USER_AGENT'} = '';
if( ($inifile ne "") && ($template ne "") ) {
    my $child_pid = fork;
    if($child_pid){
        while(! (-e $instancexml)){
            sleep 3;
        }
        print redirect(-uri=>"./view_workflow_pipeline.cgi?instance=$instancexml");    
        exit;
    } else {
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
            ## open the debugging file if needed
            open (my $debugfh, ">>$debug_file") if $debugging;
            
            chdir '/usr/local/scratch/workflow';
            use POSIX qw(setsid);
            setsid() or die "Can't start a new session: $!";
            
            setup_environment();
            
            $ENV{'PATH'} = "$ENV{'WF_ROOT'}:$ENV{'WF_ROOT'}/bin:$ENV{'WF_ROOT'}/add-ons/bin:$ENV{'PATH'}";
            
            print $debugfh $ENV{'PATH'} if $debugging;
            my $createexecstr = "$ENV{'WF_ROOT'}/CreateWorkflow -debug -i $instancexml -t $template -c $inifile --delayedbuild=true --autobuild=false > $instancexml.create.out 2>&1";
            system("$createexecstr");
            my $runexecstr = "$ENV{'WF_ROOT'}/RunWorkflow -i $instancexml > $instancexml.run.out 2>&1";
            system($runexecstr);
            close $debugfh if $debugging;
            exit;
        }
        
        exit;
    }
} else {
    my $child_pid = fork;
    if($child_pid){
        while(! (-e $instancexml)){
            sleep 3;
        }
        print redirect(-uri=>"./view_workflow_pipeline.cgi?instance=$instancexml");    
        exit;
    } else {
        close STDOUT;
        close STDERR;
        close STDIN;
        
        #  Fork again.  This helps separate the background process from
        #  the httpd process.  If we're in the original child, $gpid will
        #  hold the process id of the "grandchild", and if we're in the
        #  grandchild it will be zero.
        my $gpid = fork;
        if (! $gpid){
            ## open the debugging file if needed
            open (my $debugfh, ">>$debug_file") if $debugging;
        
            chdir '/usr/local/scratch/workflow';
            use POSIX qw(setsid);
            setsid() or die "Can't start a new session: $!";
            
            setup_environment();
            
            my $runstring = "$ENV{'WF_ROOT'}/RunWorkflow -i $instancexml -l $instancexml.log -v5 >& $instancexml.run.out";

            if ($debugging) {
                #use Data::Dumper;
                #print $debugfh Dumper(\%ENV);
                print $debugfh "preparing to run $runstring\n";
            }

            my $rc = 0xffff & system($runstring);

            printf $debugfh "system(%s) returned %#04x: $rc" if $debugging;
            if($rc == 0) {
                print $debugfh "ran with normal exit\n" if $debugging;
            } elsif ( $rc == 0xff00 ) {
                print $debugfh "command failed: $!\n" if $debugging;
            } elsif (($rc & 0xff) == 0) {
                $rc >>= 8;
                print $debugfh "ran with non-zero exit status $rc\n" if $debugging;
            } else {
                print $debugfh "ran with " if $debugging;
                if($rc & 0x80){
                    $rc &= ~0x80;
                    print $debugfh "coredump from " if $debugging;
                }
                print $debugfh "signal $rc\n" if $debugging;
            }
            
            close $debugfh if $debugging;
        }
        
        exit;
    }
}

sub setup_environment {

    umask 000;

    ## remove the apache SERVER variables from the environment
    for my $k (keys %ENV) {
        if ($k =~ /^SERVER_/ ) {
            delete $ENV{$k};
        }
    }

    ## do these need to remain, or will workflow handle them later?
    $ENV{SGE_ROOT} = '/local/n1ge';
    $ENV{SGE_CELL} = 'tigr';
    $ENV{SGE_QMASTER_PORT} = '536';
    $ENV{SGE_EXECD_PORT} = '537';
    $ENV{SGE_ARCH} = 'lx26-x86';

    $ENV{'WF_ROOT'} = "/usr/local/devel/ANNOTATION/workflow-2.2";
    #$ENV{'WF_ROOT'} = "/usr/local/devel/ANNOTATION/workflow-test-sge";
    $ENV{'WF_ROOT_INSTALL'} = "/usr/local/devel/ANNOTATION/workflow-2.2";
    #$ENV{'WF_ROOT_INSTALL'} = "/usr/local/devel/ANNOTATION/workflow-test-sge";

    $ENV{'SYBASE'} = "/usr/local/packages/sybase";
    $ENV{'PATH'} = "$ENV{'WF_ROOT'}:$ENV{'WF_ROOT'}/bin:$ENV{'WF_ROOT'}/add-ons/bin:$ENV{'PATH'}";
    $ENV{'LD_LIBRARY_PATH'} = "/usr/local/lib:/usr/local/packages/sybase/OCS/lib:/usr/local/packages/sybase/lib";
}

exit;

#use lib "/home/condor/production/request/lib";
#use TIGR::HTCRequest;


#print "hello\n";

#my $request = TIGR::HTCRequest->new(group => "antware",
#                    initialdir => "/usr/local/annotation/FIXGHI/BER2multi_dir/WORKING");
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

