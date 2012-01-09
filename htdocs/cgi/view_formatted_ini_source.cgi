#!/usr/bin/perl -w

use strict;
use CGI;
use CGI::Carp qw(fatalsToBrowser);
use File::Basename;
use Ergatis::Common;
use Ergatis::ConfigFile;
use HTML::Template;

my $q = new CGI;

print $q->header( -type => 'text/html' );

my $tmpl = HTML::Template->new( filename => 'templates/view_formatted_ini_source.tmpl',
                                die_on_bad_params => 1,
                              );

## read the ergatis config file
my $ergatis_cfg = new Ergatis::ConfigFile( -file => "ergatis.ini" );

## Pretty damn messy way of ensuring user can view non-protected files (i.e. ergatis.ini)
## when authentication is on
## TODO: Figure out a better way to allow access to non-critical configurations
my %exempt_files = ( 'ergatis.ini' => 1 );

## will be like:
## /usr/local/scratch/annotation/TTA1/workflow/runtime/wu-blastp/20724_AllGroup.niaa/component.conf.bld.ini
my $file = $q->param("file") || die "pass file";
my $pipeline_id = $q->param("pipeline_id") || undef;

## don't do it if the file doesn't end in .ini or .conf or config
if ($file !~ /\.ini$/ && $file !~ /\.conf$/ && $file !~ /\.config$/ && 
    ! $ergatis_cfg->val('display_settings', 'show_all_files')) {
    print STDERR "skipped display of $file in source viewer\n";
    quitNicely("i decline to show this type of file.");
}

validate_user_authorization($ergatis_cfg, $pipeline_id) if (! exists($exempt_files{basename($file)}));

my $sections = [];

## open the file and print it to the screen.
open (my $ifh, "<$file") || quitNicely("couldn't open file $file");

my $ini = new Ergatis::ConfigFile( -file => $file );

for my $section ( $ini->Sections ) {
    my $parameters = [];
    
    for my $parameter ( $ini->Parameters($section) ) {
        my $value = $ini->val($section, $parameter);
        my $url;
        
        ## get any comment and format
        my $comment = $ini->get_comment_html($section, $parameter);
        
        ## look for any linkable xml
        if ( $value =~ m^(?<!\$\;)(/[/a-z0-9_\-.]+\.(?:xml|instance|bsml))\s*$^i ) {
            $url = $1;
            $value =~ s|$url|<a href="./view_formatted_xml_source.cgi?file=$url&pipeline_id=$pipeline_id">$url</a>|;
        }

        ## look for any linkable ini
        if ( $value =~ m^(?<!\$\;)(/[/a-z0-9_\-.]+\.(?:ini|config|conf))\s*$^i ) {
            $url = $1;
            $value =~ s|$url|<a href="./view_formatted_ini_source.cgi?file=$url&pipeline_id=$pipeline_id">$url</a>|;
        }

        ## look for any linkable lists
        if ( $value =~ m^(?<!\$\;)(/[/a-z0-9_\-.]+\.list)\s*$^i ) {
            $url = $1;
            $value =~ s|$url|<a href="./view_raw_source.cgi?file=$url&pipeline_id=$pipeline_id">$url</a>|;
        }
        
        push @$parameters, { parameter => $parameter, value => $value, comment => $comment };
    }
    
    push @$sections, { section => $section, parameters => $parameters };
}

$tmpl->param( FILE                => $file );
$tmpl->param( SECTIONS            => $sections );

$tmpl->param( QUICK_LINKS         => &get_quick_links($ergatis_cfg) );
$tmpl->param( SUBMENU_LINKS       => [
                                        { label => 'view unformatted version', is_last => 1, url => "./view_raw_source.cgi?file=$file&pipeline_id=$pipeline_id" },
                                     ] );

print $tmpl->output;


sub quitNicely {
    my $msg = shift;
    
    print $msg;
    exit;
}
