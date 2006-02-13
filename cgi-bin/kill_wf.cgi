#!/usr/local/bin/perl

use strict;

use XML::Twig;

use CGI qw(:standard);

use File::Find;
use File::Basename;
use Ergatis::ConfigFile;
use Proc::ProcessTable;

my $instancexml = param('instancexml');

my $ergatis_cfg = new Ergatis::ConfigFile( -file => "ergatis.ini" );

if(! -d $ergatis_cfg->val('paths','workflow_run_dir')){
    die "Invalid workflow_run_dir in ergatis.ini : " . $ergatis_cfg->val('paths','workflow_run_dir');
}

my $rundir = $ergatis_cfg->val('paths','workflow_run_dir');
my $repository_root;
my $project;
my $pipelineid;

if ( $instancexml =~ m|(.+/(.+?))/Workflow/pipeline/(\d+)/| ) {
    $repository_root = $1;
    $project = $2;
    $pipelineid = $3;
} else {
    die "failed to extract a repository_root from $instancexml.  expected a Workflow subdirectory somewhere."
}

my ($pid,$hostname,$execuser) = &parselockfile("$repository_root/workflow/lock_files/pid.$pipelineid");

my $t = new Proc::ProcessTable;
my $parentpids = {};
$parentpids->{$pid} = 0;
my @childtokill;
push @childtokill,$pid;
foreach my $p ( @{$t->table} ){
    if(exists $parentpids->{$p->{ppid}}){
	$parentpids->{$p->{pid}} = $p->{ppid};
	if($p->{ppid} eq $childtokill[$#childtokill]){
	    push @childtokill,$p->{pid};
	}
    }
}
my $count;

if($childtokill[3] > 1){
    $count = kill 15, $childtokill[3];
}

if ($count != 1 || $childtokill[3] eq '') {
    print header();
    print "<html><body>";
    print "Attempting to kill process $execuser\@$hostname:$pid. Error signalling process with pid $pid; manual intervention may be required<br><a href='./view_workflow_pipeline.cgi?instance=$instancexml'>[view workflow]</a><br><i>Debug information follows</i><br>";
    print "Detected workflow processes: ",join(',',@childtokill),"<br>";
    print "<hr>";
    foreach my $p ( @{$t->table} ){
	print "Pid: $p->{pid} parent pid: $p->{ppid} cmdline: $p->{cmndline}<br>";
    }

    print "</body></html";
    exit(1);
} else {
    print redirect(-uri=>"./view_workflow_pipeline.cgi?instance=$instancexml");
    exit(0);
}
print header();
print "<html><body>Can't find process for instance $instancexml.  Manual intervention required</body></html>";
exit;

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
	return ($pid,$hostname,$getpwuid);
    }
    return undef;
}

