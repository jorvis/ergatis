#!/usr/local/bin/perl
use strict;

use CGI qw(:standard);

use Tree::DAG_Node;
use File::Find;
use File::Basename;
use Ergatis::Validator;
use Ergatis::ConfigFile;

use HTML::Template;
use Sys::Hostname;

## toggles debugging messages
my $debugging = 0;
my $debug_file = '/tmp/ergatis.debug';

my $instancexml = param('instancexml');
my $inifile = param('inifile');
my $template = param('xmltemplate');
my $validate = param('validate');
my $pipelineid = param('pipelineid');

$validate = 1 if (! defined $validate);

my $log =  $instancexml.".log";
my $out = $instancexml.".out";

my $ergatis_cfg = new Ergatis::ConfigFile( -file => "ergatis.ini" );

if(! -d $ergatis_cfg->val('paths','workflow_run_dir')){
    die "Invalid workflow_run_dir in ergatis.ini : " . $ergatis_cfg->val('paths','workflow_run_dir');
}
if(! -d $ergatis_cfg->val('paths','workflow_root')){
    die "Invalid workflow_root in ergatis.ini : " . $ergatis_cfg->val('paths','workflow_root');
}
if(! -e $ergatis_cfg->val('paths','workflow_log4j')){
    die "Invalid workflow_log4j in ergatis.ini : " . $ergatis_cfg->val('paths','workflow_log4j');
}

## extract the repository root from the instancexml path
my $repository_root;
if ( $instancexml =~ m|^(.+)/Workflow| ) {
    $repository_root = $1;
} else {
    die "expected a Workflow dir in the instancexml path";
}

my $lockfile = "$repository_root/workflow/lock_files/pid.$pipelineid";

if(!$pipelineid){
    die "Required parameter --pipelineid missing\n";
}

if ( $validate ) {

    my $validator = new Ergatis::Validator;
    my $warnings = $validator->check_pipeline_template( pipeline => $template );

    if ( scalar @$warnings ) {
        print "Content-type: text/html\n\n";

        my $tmpl = HTML::Template->new( filename => 'templates/warning_intermediate.tmpl',
                                        die_on_bad_params => 1,
                                      );
        
        my $warning_hash = [];
        for my $warning ( @$warnings ) {
            push @$warning_hash, { text => $warning };
        }

        $tmpl->param( MESSAGE => 'The following are potential problems detected with your pipeline.  Use the back button ' .
                                 'on your browser if you need to go back and fix these, or click the \'continue anyway\' ' .
                                 'button to ignore them.' );
        $tmpl->param( WARNING_MESSAGES => $warning_hash );
        $tmpl->param( CONTINUE_ACTION  => "./run_wf.cgi?xmltemplate=$template&amp;inifile=$inifile&amp;instancexml=$instancexml&amp;validate=0&amp;pipelineid=$pipelineid" );
        $tmpl->param( RETRY_ACTION => "./run_wf.cgi?xmltemplate=$template&amp;inifile=$inifile&amp;instancexml=$instancexml&amp;validate=1&amp;pipelineid=$pipelineid" );

        print $tmpl->output;
        exit;
    }
}

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
            
            chdir $ergatis_cfg->val('paths','workflow_run_dir');
            use POSIX qw(setsid);
            setsid() or die "Can't start a new session: $!";
            
            setup_environment();
            
            $ENV{'PATH'} = "$ENV{'WF_ROOT'}:$ENV{'WF_ROOT'}/bin:$ENV{'WF_ROOT'}/add-ons/bin:$ENV{'PATH'}";
            
            print $debugfh $ENV{'PATH'} if $debugging;
            my $createexecstr = "$ENV{'WF_ROOT'}/CreateWorkflow -debug -i $instancexml -t $template -c $inifile --delayedbuild=true --autobuild=false > $instancexml.create.out 2>&1";
            system("$createexecstr");

            my $runexecstr = "$ENV{'WF_ROOT'}/RunWorkflow -i $instancexml --logconf=".$ergatis_cfg->val('paths','workflow_log4j')." > $instancexml.run.out 2>&1";

        &writelockfile($lockfile);
            
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
        
            chdir $ergatis_cfg->val('paths','workflow_run_dir');
            use POSIX qw(setsid);
            setsid() or die "Can't start a new session: $!";
            
            setup_environment();

            my $runstring = "$ENV{'WF_ROOT'}/RunWorkflow -i $instancexml -l $instancexml.log --logconf=".$ergatis_cfg->val('paths','workflow_log4j')." >& $instancexml.run.out";
        #
        #------------------------------------------------------------------------------------------------------------

            if ($debugging) {
                #use Data::Dumper;
                #print $debugfh Dumper(\%ENV);
                print $debugfh "preparing to run $runstring\n";
            }
        
        &writelockfile($lockfile);

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

    $ENV{'WF_ROOT'} = $ergatis_cfg->val('paths','workflow_root');
    $ENV{'WF_ROOT_INSTALL'} = $ergatis_cfg->val('paths','workflow_root');

    $ENV{'SYBASE'} = "/usr/local/packages/sybase";
    $ENV{'PATH'} = "$ENV{'WF_ROOT'}:$ENV{'WF_ROOT'}/bin:$ENV{'WF_ROOT'}/add-ons/bin:$ENV{'PATH'}";
    $ENV{'LD_LIBRARY_PATH'} = "";
}

sub writelockfile{
    my($file) = @_;
    my $retrycount=0;
    if(-e $file){
    my($pid,$hostname,$getpwuid,$retries) = &parselockfile($file);
    $retries=0 if(!$retries);
    $retrycount = $retries++;
    }
    open FILE,"+>$file" or die "Can't open lock file $file for writing :$!";
    print FILE $$,"\n";
    print FILE hostname(),"\n";
    my $user = getpwuid($<);
    print FILE $user,"\n";
    print FILE $retrycount,"\n";
    close FILE;
}

sub parselockfile{
    my($file) = @_;
    if(-e $file){
    open FILE, "$file" or die "Can't open lock file $file";
    my(@elts) = <FILE>;
    close FILE;
    chomp(@elts);
    my $pid = $elts[0];
    my $hostname = $elts[1];
    my $getpwuid = $elts[2];
    my $retries = $elts[3];
    return ($pid,$hostname,$getpwuid,$retries);
    }
    return undef;
}


exit;

