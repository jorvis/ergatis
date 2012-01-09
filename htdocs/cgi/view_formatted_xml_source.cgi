#!/usr/bin/perl -w

use strict;
use CGI;
use CGI::Carp qw(fatalsToBrowser);
use Ergatis::Common;
use Ergatis::ConfigFile;
use File::Basename;
use HTML::Template;

my $q = new CGI;

print $q->header( -type => 'text/html' );

my $tmpl = HTML::Template->new( filename => 'templates/view_formatted_xml_source.tmpl',
                                die_on_bad_params => 1,
                              );

## read the ergatis config file
my $ergatis_cfg = new Ergatis::ConfigFile( -file => "ergatis.ini" );

## will be like:
## /usr/local/scratch/annotation/TGA1/workflow/runtime/split_fasta/29134_test2/pipeline.xml
my $file = $q->param("file") || die "pass file";
my $pipeline_id = $q->param("pipeline_id") || undef;

## the file may have been compressed
if ( ! -e $file && -e "$file.gz" ) {
    $file .= '.gz';
}

## don't do it if the file doesn't end in .xml or .instance
if ($file !~ /\.xml$/ && 
    $file !~ /\.instance$/ && 
    $file !~ /\.bsml$/ &&
    $file !~ /\.gz$/ &&
    ! $ergatis_cfg->val('display_settings', 'show_all_files')) {
    print STDERR "skipped display of $file in source viewer\n";
    quitNicely("i decline to show this type of file.");
}

## Grab the repository root from our file path
my $repository_root = "unknown";
my $project;

## If we aren't passed in an ID attempt to grab it from the path
if ( !($pipeline_id) && $file =~ m|(.+/(.+?))/workflow/runtime/.+?/([A-Z0-9]+)(_.+)?/| ) { 
    $repository_root = $1; 
    $project = $2;
    $pipeline_id = $3;
} 

## If per-account pipeline security is enabled we will want to ensure that the user currently logged in
## has access to this pipeline.
validate_user_authorization($ergatis_cfg, $pipeline_id);

my $progress_image_width = 500;

## open the file and print it to the screen.
my $ifh;
if ($file =~ /\.gz$/) {
    open($ifh, "<:gzip", $file) || quitNicely("couldn't open file $file");
} else {
    open($ifh, "<$file") || quitNicely("couldn't open file $file");
}


## something to remember what states we find
my %states;
my $overall_state = 0;
my $found_states = 0;
my $has_multiple_states = 0;
my $within_status_box = 0;
my $command_count = 0;
#my @xmlfiles;
my %xmlfiles;
my $display_source = '';

while (my $line = readline $ifh) {
    my ($tag, $url, $word);

    ## have we found the first state?
    if ( (! $overall_state) && $line =~ m|<state>(.+)</state>|) {
        $overall_state = $1;
    }

    ## are we dealing with status box?
    if ($line =~ /<status>/) {
        $within_status_box = 1;
    } elsif ($line =~ /<\/status>/) {
        $within_status_box = 0;
        $found_states = 1;
    } elsif ($within_status_box && (! $found_states) && $line =~ /\<(.+?)\>(\d+)/) {
        my $state = $1;
        $state = 'error' if $state eq 'errors';  ## workflow doesn't always match states exactly.
    
        if ($2 && $state ne 'total') {
            $states{$state} += $2;
            $command_count += $2;
        }
    }
    
    ## colorize any tags
    my $newline = $line;
    while ( $line =~ m|\<(/*\w+)|g ) {
        $tag = $1;
        $newline =~ s|$tag|STARTTAGSPAN${tag}ENDTAGSPAN|;
    }

    $line = $newline;
    
    ## look for comments
    $line =~ s/<!--/STARTCOMMENTSPAN<!--/g;
    $line =~ s|-->|-->ENDCOMMENTSPAN|g;

    $line =~ s/\</\&lt\;/g;
    $line =~ s/\>/\&gt\;/g;

    ## add the spans properly
    $line =~ s|STARTTAGSPAN|<span class="tag">|g;
    $line =~ s|STARTCOMMENTSPAN|<span class="comment">|g;
    $line =~ s!(ENDTAGSPAN|ENDCOMMENTSPAN)!</span>!g;

    ## look for any linkable xml
    if ( $line =~ m^(?<!\$\;)(/[/a-z0-9_\-.]+\.(?:xml|instance|bsml))^i ) {
        $url = $1;
        $line =~ s|$url|<a href="./view_formatted_xml_source.cgi?pipeline_id=$pipeline_id&file=$url">$url</a>|;
        $xmlfiles{$url}++;
    }
    
    ## look for any linkable ini
    if ( $line =~ m^(?<!\$\;)(/[/a-z0-9_\-.]+\.(?:ini|config|conf))^i ) {
        $url = $1;
        $line =~ s|$url|<a href="./view_formatted_ini_source.cgi?pipeline_id=$pipeline_id&file=$url">$url</a>|;
    }

    ## look for any linkable log/stderr/sdtout
    if ( $line =~ m^(?<!\$\;)(/[/a-z0-9_\-.]+\.(?:log|stderr|stdout))^i ) {
        if(-z $url){
            $url = $1;
            $line =~ s|$url|<a href="./view_formatted_log_source.cgi?pipeline_id=$pipeline_id&file=$url">$url</a>|;
        }
    }

    ## look for any linkable lists
    if ( $line =~ m^(?<!\$\;)(/[/a-z0-9_\-.]+\.list)\s*$^i ) {
        $url = $1;
        $line =~ s|$url|<a href="./view_raw_source.cgi?pipeline_id=$pipeline_id&file=$url">$url</a>|;
    }

    ##match any other files
    if ( $line =~ m|(?<!\$\;)(/[/a-z0-9_\-\.]+)\&|i ) {
        $url = $1;
        if(-f $url){
            $line =~ s|$url|<a href="./view_raw_source.cgi?pipeline_id=$pipeline_id&file=$url">$url</a>|;
        }
    }

    ## look for any execution hosts
    if ( $line =~ m|executionHost</span>\&gt\;(.+?)&lt;|i ) {
        my $ehost = $1; 
        my $hostsrvstr = join(',',split(/\./,$ehost));
        $line =~ s|$ehost|<a href="http://intranet.tigr.org/cgi-bin/sysadmin/hobbit/bb-hostsvc.sh?HOSTSVC=$hostsrvstr.cpu&amp;IP=0.0.0.0&amp;DISPLAYNAME=$ehost">$ehost</a>|;
    }
    
    $display_source .= $line;
}

my $state_elements = [];
my $linked_files = [];
my $list_limit = 10;
my $unshown_file_count = 0;

if (scalar keys %states) {

    ## build the line that lists each status and its count
    ## at the same time we can calculate the width of each state for the progress bar
    my $state_count = scalar keys %states;
    my $width_used = 0;  
    my $states_handled = 0;
   
    for my $status (sort keys %states) {
        ## each status gives a percentage of the total command_count
        my $width = int( ($states{$status} / $command_count) * $progress_image_width);
        $width_used += $width;
        $states_handled++;
        
        ## if this is the last state and there are unused pixels in the bar, just tag
        ## the unused ones onto this color so we don't have a gap
        if ( $states_handled == $state_count && $width_used < $progress_image_width ) {
            $width += $progress_image_width - $width_used;
        }

        push @$state_elements, { state => $status, count => $states{$status}, width => $width };
    }

    if ($state_count > 1) {
        $has_multiple_states = 1;
    }

    my $file_counter = 0;
    for (keys %xmlfiles) {
        if ( $file_counter == $list_limit ) {
            $unshown_file_count = scalar( keys %xmlfiles ) - $file_counter;
            last;
        }
        
        my @stats = stat $_; 
        
        push @$linked_files, {
                                url => "./view_formatted_xml_source.cgi?pipeline_id=$pipeline_id&file=$_",
                                label => basename($_),
                                size => sprintf("%.1f", $stats[7]/1024) . ' kb',
                             };
    }
}

## If we were able to parse a project out of the XML file we are viewing add a link back to the projects
## pipeline list
my @submenu_links = ();
if ($repository_root ne "unknown") {
        push(@submenu_links, { label => 'pipeline list', is_last => 0, url => "./pipeline_list.cgi?repository_root=$repository_root" })
}
push(@submenu_links, { label => 'view unformatted version', is_last => 1, url => "./view_raw_source.cgi?file=$file&pipeline_id=$pipeline_id" });

$tmpl->param( FILE                => $file );
$tmpl->param( DISPLAY_SOURCE      => $display_source );
$tmpl->param( STATE_ELEMENTS      => $state_elements );
$tmpl->param( LINKED_FILES        => $linked_files );
$tmpl->param( OVERALL_STATE       => $overall_state );
$tmpl->param( HAS_MULTIPLE_STATES => $has_multiple_states );
$tmpl->param( UNSHOWN_FILE_COUNT  => $unshown_file_count );

$tmpl->param( QUICK_LINKS         => &get_quick_links($ergatis_cfg) );
$tmpl->param( SUBMENU_LINKS       => \@submenu_links);

print $tmpl->output;


sub quitNicely {
    my $msg = shift;
    
    print $msg;
    exit;
}
