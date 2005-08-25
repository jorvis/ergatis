#!/usr/local/bin/perl -w

use strict;
use CGI;
use CGI::Carp qw(fatalsToBrowser);
use Date::Manip;
use File::stat;
use POSIX;
use XML::Twig;

my $q = new CGI;
my $repository_root;

print $q->header( -type => 'text/html' );

## will be like:
## /usr/local/scratch/annotation/TGA1/Workflow/split_fasta/29134_test2/pipeline.xml
my $pipeline = $q->param("pipeline") || die "pass pipeline";

## get the project label out of the pipeline.  just before 'Workflow'
my $project = '?';
if ( $pipeline =~ m|.+/(.+?)/Workflow| ) {
    $project = $1;
}


my $twig = new XML::Twig;

if ($pipeline =~ /(.+)\/Workflow/ ) {
    $repository_root = $1;
}

$twig->parsefile($pipeline);

my $commandSetRoot = $twig->root;
my $commandSet = $commandSetRoot->first_child('commandSet');

## pull desired info out of the root commmandSet
## pull: project space usage
my ($starttime, $endtime, $lastmodtime, $state, $runtime) = ('n/a', 'n/a', '', 'unknown', 'n/a');
my ($starttimeobj, $endtimeobj);

if ($commandSet->first_child('startTime') ) {
    $starttimeobj = ParseDate($commandSet->first_child('startTime')->text());
    $starttime = UnixDate($starttimeobj, "%c");
}

if ($commandSet->first_child('endTime') ) {
    $endtimeobj = ParseDate($commandSet->first_child('endTime')->text());
    $endtime = UnixDate($endtimeobj, "%c");
}

if ( $commandSet->first_child('state') ) {
    $state  = $commandSet->first_child('state')->text();
}

## we can calculate runtime only if start and end time are known, or if start is known and state is running
if ($starttimeobj) {
    if ($endtimeobj) {
        $runtime = strftime( "%H hr %M min %S sec", reverse split(/:/, DateCalc($starttimeobj, $endtimeobj)) );
    } elsif ($state eq 'running') {
        $runtime = strftime( "%H hr %M min %S sec", reverse split(/:/, DateCalc("now", $starttimeobj)) ) . ' ...';
    }
}

my $filestat = stat($pipeline);
my $user = getpwuid($filestat->uid);
$lastmodtime = time - $filestat->mtime;
$lastmodtime = strftime( "%H hr %M min %S sec", reverse split(/:/, DateCalc("today", ParseDate(DateCalc("now", "- ${lastmodtime} seconds")) ) ));

## quota information (only works if in /usr/local/annotation/SOMETHING)
my $quotastring = 'quota information currently disabled';
#my $quotastring = '';
#if ($pipeline =~ m|^(/usr/local/annotation/.+?/)|) {
#    $quotastring = `getquota -N $1`;
#    if ($quotastring =~ /(\d+)\s+(\d+)/) {
#        my ($limit, $used) = ($1, $2);
#        $quotastring = sprintf("%.1f", ($used/$limit) * 100) . "\% ($used KB of $limit KB used)";
#    } else {
#        $quotastring = 'error parsing quota information';
#    }
#} else {
#    $quotastring = 'unavailable (outside of project area)';
#}

print <<PipelineSummarY;
    <div id='pipeline'>$pipeline</div>
    <div class='pipelinestat' id='pipelinestart'><strong>start:</strong> $starttime</div>
    <div class='pipelinestat' id='pipelineend'><strong>end:</strong> $endtime</div>
    <div class='pipelinestat' id='pipelinelastmod'><strong>last mod:</strong> $lastmodtime</div><br>
    <div class='pipelinestat'><strong>state:</strong> <span id='pipelinestate'>$state</span></div>
    <div class='pipelinestat' id='pipelineuser'><strong>user:</strong> $user</div>
    <div class='pipelinestat' id='pipelineruntime'><strong>runtime:</strong> $runtime</div><br>
    <div class='pipelinestat'><strong>project:</strong> <span id='projectid'>$project</span></div>
    <div class='pipelinestat' id='projectquota'><strong>quota:</strong> $quotastring</div>
    <div class='timer' id='pipeline_timer_label'></div>
    <div id='pipelinecommands'>
        [<a href='./new_pipeline.cgi?&root=$repository_root/Workflow/pipeline'>new</a>] 
        [<a href='./run_wf.cgi?instancexml=$pipeline'>rerun</a>] 
        [<a href='./show_pipeline.cgi?xmltemplate=$pipeline&edit=1'>edit</a>] 
        [<a href='./kill_wf.cgi?instancexml=$pipeline'>kill</a>] 
        [<a href='http://htcmaster.tigr.org/antware/condor-status/index.cgi' target='_blank'>condor status</a>] 
        [<a href='http://intranet.tigr.org/grid/cgi-bin/sgestatus.cgi' target='_blank'>SGE status</a>] 
    </div>
PipelineSummarY







