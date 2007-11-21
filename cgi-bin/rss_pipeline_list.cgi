#!/usr/bin/perl -w

use strict;
use CGI;
use CGI::Carp qw(fatalsToBrowser);
use Ergatis::Common;
use Ergatis::ConfigFile;
use File::Basename;
use Storable;
use XML::RSS;

my $q = new CGI;

## serve mozilla tex/xml, others application/rss+xml
##  strategy described here: http://www.petefreitag.com/item/381.cfm
if ($q->user_agent() =~ /mozilla/i) {
    print $q->header( -type => 'text/xml' );
} else {
    print $q->header( -type => 'application/rss+xml' );  ## for rss 2.0
}

## read the ergatis config file
my $ergatis_cfg = new Ergatis::ConfigFile( -file => "ergatis.ini" );

my $temp_space = $ergatis_cfg->val( 'paths', 'temp_space' ) || die "temp_space not defined in ergatis.ini file";

## an MD5 is used on the project list so that many installations can
## share the same storable object.
my $cfg_md5 = $ergatis_cfg->project_list_md5();

my $running_pipelines = [];
my $active_pipelines = [];

## is there a stored version of the running pipeline list?
my $running_dump_file = "$temp_space/$cfg_md5.ergatis.running.dump";
my $active_dump_file  = "$temp_space/$cfg_md5.ergatis.active.dump";

if ( -e $running_dump_file && -e $active_dump_file) {
    ## pull the lists from the disk.
    $running_pipelines = retrieve $running_dump_file;
    $active_pipelines  = retrieve $active_dump_file;
}

my $base_url = $q->url();

my $rss = new XML::RSS (version => '2.0');

$rss->channel(  title       => 'ergatis',
                link        => $base_url,
                language    => 'en',
                description => 'ergatis active pipelines',
#                pubDate     => '',  # like 'Thu, 23 Aug 1999 07:00:00 GMT'
#                lastBuildDate => '',
             );

for my $pipeline ( @$running_pipelines, @$active_pipelines ) {
    my $url = "$base_url/$pipeline->{pipeline_url}";
    
    ## take of the current script part of the path
    $url =~ s|rss_pipeline_list\.cgi/\./||;
    
    $rss->add_item(
        title       => "$pipeline->{pipeline_id} | $pipeline->{label} | $pipeline->{state}",
        link        => $url,
        description => "pipeline $pipeline->{pipeline_id} in project $pipeline->{label} is currently $pipeline->{state}",
#        pubDate     => 'Thu, 6 Apr 2006 04:09:36 GMT',
    );
}

print $rss->as_string();
