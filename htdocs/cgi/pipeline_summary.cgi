#!/usr/bin/perl -w

use strict;
use CGI;
use CGI::Carp qw(fatalsToBrowser);
use Date::Manip;
use Ergatis::Monitor;
use HTML::Template;
use XML::Twig;

my $q = new CGI;
my $repository_root;

print $q->header( -type => 'text/html' );

my $tmpl = HTML::Template->new( filename => 'templates/pipeline_summary.tmpl',
                                die_on_bad_params => 1,
                              );

my $pipeline = $q->param("pipeline") || die "pass pipeline";
my $project;
my $pipelineid;
if ( $pipeline =~ m|(.+/(.+?))/workflow/runtime/pipeline/([A-Z0-9]+)/| ) {
    $repository_root = $1;
    $project = $2;
    $pipelineid = $3;
} else {
    die "failed to extract a repository_root from $pipeline.  expected a workflow/runtime subdirectory somewhere."
}

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
    
    ## if the state is 'incomplete', check for a token file that indicates
    #   that this pipeline was submitted to a job manager.  this allows us to
    #   show a 'pending' state of the parent pipeline before the XML is parsed.
    if ( $state eq 'incomplete' && -e "$pipeline.submitted" ) {
        $state = 'pending';
    }
}

my ($starttime, $endtime, $runtime) = &time_info( $commandSet );

## endtime can be null
$endtime = '?' unless $endtime;

my $filestat = stat($pipeline);
my $user = getpwuid($filestat->uid);
my $lastmodtime = time - $filestat->mtime;
$lastmodtime = 0 if( $lastmodtime < 0 );
$lastmodtime = strftime( "%H hr %M min %S sec", reverse split(/:/, DateCalc("today", ParseDate(DateCalc("now", "- ${lastmodtime} seconds")) ) ));

## quota information (only works if in /usr/local/annotation/SOMETHING)
my $quotastring = &quota_string($repository_root);

$tmpl->param( PIPELINE_FILE       => $pipeline );
$tmpl->param( START_TIME          => $starttime );
$tmpl->param( END_TIME            => $endtime );
$tmpl->param( LAST_MOD_TIME       => $lastmodtime );
$tmpl->param( PIPELINE_STATE      => $state );
$tmpl->param( USER                => $user );
$tmpl->param( RUNTIME             => $runtime );
$tmpl->param( PROJECT             => $project );
$tmpl->param( QUOTA_STRING        => $quotastring );
$tmpl->param( PIPELINE_ID         => $pipelineid );

print $tmpl->output;
