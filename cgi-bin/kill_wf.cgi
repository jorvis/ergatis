#!/usr/local/bin/perl

use strict;

use XML::Twig;

use CGI qw(:standard);

use File::Find;
use File::Basename;

my $instancexml = param('instancexml');

my $rundir = '/usr/local/scratch/workflow';
opendir(DIR,$rundir) or die "Can't open directory $rundir";
while(defined (my $file = readdir(DIR))){
    if($file =~ /^.workflow.*\d+$/){
	my($pid,$invocation) = &parse_run_file("$rundir/$file");
	if($invocation =~ /$instancexml/){
	    my $count = kill 15, $pid;
	    if ($count != 1) {
		print header();
		print "<html><body>";
		print "Error signalling process with pid $pid; manual intervention may be required<br><a href='./view_workflow_pipeline.cgi?xml_input=$instancexml'>[view workflow]</a>";
		print "</body></html";
		exit(1);
	    } else {
		print redirect(-uri=>"./view_workflow_pipeline.cgi?xml_input=$instancexml");
		exit(0);
	    }
	}
    }
}
closedir DIR;
print header();
print "<html><body>Can't find process for instance $instancexml.  Manual intervention required</body></html>";
exit;


sub parse_run_file{
    my ($file) = @_;
    open FILE, $file or die "Can't open file $file\n";
    my $pid;
    my $invocation;
    while(my $line=<FILE>){
	if($line =~ /^pid/){
	    ($pid) = ($line =~ /pid\s*=\s*(\d+)/);
	}
	elsif($line =~ /^invocation/){
	    ($invocation) = ($line =~ /invocation\s*=\s*(.+)/);
	}
    }
    close FILE;
    return ($pid,$invocation);
}
