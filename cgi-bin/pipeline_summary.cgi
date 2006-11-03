#!/usr/local/bin/perl -w

use strict;
use CGI;
use CGI::Carp qw(fatalsToBrowser);
use Date::Manip;
use File::stat;
use Monitor;
use POSIX;
use XML::Twig;

my $q = new CGI;
my $repository_root;

print $q->header( -type => 'text/html' );

## will be like:
## /usr/local/scratch/annotation/TGA1/Workflow/split_fasta/29134_test2/pipeline.xml

my $pipeline = $q->param("pipeline") || die "pass pipeline";
my $project;
my $pipelineid;
if ( $pipeline =~ m|(.+/(.+?))/workflow/runtime/pipeline/(\d+)/| ) {
    $repository_root = $1;
    $project = $2;
    $pipelineid = $3;
} else {
    die "failed to extract a repository_root from $pipeline.  expected a workflow/runtime subdirectory somewhere."
}

my ($pid,$hostname,$execuser,$retries) = &parselockfile("$repository_root/workflow/lock_files/pid.$pipelineid");
$pid = "" if(!$pid);
$execuser = "unknown" if(!$execuser);
$hostname = "unknown" if(!$hostname);

my $twig = new XML::Twig;

if ($pipeline =~ /(.+)\/workflow\/runtime/ ) {
    $repository_root = $1;
}

my $pipeline_fh;
if ($pipeline =~ /\.gz/) {
    open($pipeline_fh, "<:gzip", "$pipeline") || die "can't read $pipeline: $!"; 
} elsif ( ! -e $pipeline && -e "$pipeline.gz" ) {
    open($pipeline_fh, "<:gzip", "$pipeline.gz") || die "can't read $pipeline: $!"; 
} else {
    open($pipeline_fh, "<$pipeline") || die "can't read $pipeline: $!";       
}

$twig->parse($pipeline_fh);

my $commandSetRoot = $twig->root;
my $commandSet = $commandSetRoot->first_child('commandSet');

## pull desired info out of the root commmandSet
## pull: project space usage

my $state = 'unknown';
if ( $commandSet->first_child('state') ) {
    $state  = $commandSet->first_child('state')->text();
}

my ($starttime, $endtime, $runtime) = &time_info( $commandSet );

## endtime can be null
$endtime = '?' unless $endtime;

my $filestat = stat($pipeline);
my $user = getpwuid($filestat->uid);
my $lastmodtime = time - $filestat->mtime;
$lastmodtime = strftime( "%H hr %M min %S sec", reverse split(/:/, DateCalc("today", ParseDate(DateCalc("now", "- ${lastmodtime} seconds")) ) ));

## quota information (only works if in /usr/local/annotation/SOMETHING)
my $quotastring = &quota_string($repository_root);

print <<PipelineSummarY;
    <div id='pipeline'>$pipeline</div>
        <div class='pipelinestat' id='pipelinestart'><strong>start:</strong> $starttime</div>
        <div class='pipelinestat' id='pipelineend'><strong>end:</strong> $endtime</div>
        <div class='pipelinestat' id='pipelinelastmod'><strong>last mod:</strong> $lastmodtime</div><br>
        <div class='pipelinestat'><strong>state:</strong> <span id='pipelinestate'>$state</span></div>
        <div class='pipelinestat' id='pipelineuser'><strong>user:</strong> $user</div>
        <div class='pipelinestat' id='pipelineruntime'><strong>runtime:</strong> $runtime</div>
        <div class='pipelinestat' id='pipelineretry'><strong>retries:</strong> $retries</div>
        <div class='pipelinestat' id='pipelineexec'><strong>exec host:</strong> <a href='http:$hostname:8080/ergatis/view_pipeline.cgi?instance=$pipeline'>$execuser\@$hostname</a>:$pid</div><br>
        <div class='pipelinestat'><strong>project:</strong> <span id='projectid'>$project</span></div>
        <div class='pipelinestat' id='projectquota'><strong>quota:</strong> $quotastring</div>
        <div class='pipelinestat' id='pipelineid'><strong>pipeline id:</strong> $pipelineid</div>
    <div class='timer' id='pipeline_timer_label'></div>
    <div id='pipelinecommands'>
        <a href='./pipeline_list.cgi?repository_root=$repository_root'><img class='navbutton' src='/ergatis/button_blue_pipeline_list.png' alt='pipeline list' title='pipeline list'></a>
        <a href='./new_pipeline.cgi?repository_root=$repository_root'><img class='navbutton' src='/ergatis/button_blue_new.png' alt='new' title='new'></a>
        <a href='./run_wf.cgi?instancexml=$pipeline&validate=0&pipelineid=$pipelineid'><img class='navbutton' src='/ergatis/button_blue_rerun.png' alt='rerun' title='rerun'></a>  
        <a href='./show_pipeline.cgi?xmltemplate=$pipeline&edit=1'><img class='navbutton' src='/ergatis/button_blue_edit.png' alt='edit' title='edit'></a>
        <a href='./kill_wf.cgi?instancexml=$pipeline'><img class='navbutton' src='/ergatis/button_blue_kill.png' alt='kill' title='kill'></a>
        <a href='http://htc.tigr.org/antware/cgi-bin/sgestatus.cgi'><img class='navbutton' src='/ergatis/button_blue_grid_info.png' alt='grid info' title='grid info'></a>
        <a href='/cgi-bin/ergatis/view_formatted_xml_source.cgi?file=$pipeline'><img class='navbutton' src='/ergatis/button_blue_xml.png' alt='View XML' title='View XML'></a>
    </div>
PipelineSummarY







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
