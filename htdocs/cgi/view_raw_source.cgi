#!/usr/bin/perl -w

use strict;
use CGI;
use CGI::Carp qw(fatalsToBrowser);
use Ergatis::ConfigFile;
use Ergatis::Common;

my $q = new CGI;

print $q->header( -type => 'text/plain' );

## will be like:
my $file = $q->param("file") || die "pass file";


## the file may have been compressed
if ( ! -e $file && -e "$file.gz" ) {
    $file .= '.gz';
}

## read the ergatis config file
my $ergatis_cfg = new Ergatis::ConfigFile( -file => "ergatis.ini" );

## only print certain file types (for minor security)
if ($file !~ /\.xml$/ && 
    $file !~ /\.instance$/ && 
    $file !~ /\.bsml$/ &&
    $file !~ /\.gz$/ && 
    $file !~ /\.ini/ &&
    $file !~ /\.log/ &&
    $file !~ /\.config/ &&
    $file !~ /\.list/ &&
    $file !~ /\.stderr/ &&
    $file !~ /\.stdout/ && 
    ! $ergatis_cfg->val('display_settings', 'show_all_files')) {
    print STDERR "skipped display of $file in source viewer\n";
    quitNicely("i decline to show this type of file.");
}

## If our file is a component of a specific pipeline we will want to restrict 
## access to it if account pipeline security is enabled
if ( $file =~ m|(.+/(.+?))/workflow/runtime/(.+?)/([A-Z0-9]+)_(.+?)/| ) {
    my $repository_root = $1;
    my $project = $2;
    my $component_name = $3;
    my $pipeline_id = $4;
    my $output_token = $5;

    ## If per-account pipeline security is enabled we will want to ensure that the user currently logged in
    ## has access to this pipeline.
    my $auth = validate_user_authorization($ergatis_cfg, $project, $repository_root, $pipeline_id, 1);
    quitNicely("User does not have authorization to view this resource") if (! $auth);
}


## open the file and print it to the screen.
my $ifh;
if ($file =~ /\.gz$/) {
    open($ifh, "<:gzip", $file) || quitNicely("couldn't open compressed file $file: $!");
} else {
    open($ifh, "<$file") || quitNicely("couldn't open file $file");
}

while (my $line = readline $ifh) {
    ## just print it    
    print $line;
}

sub quitNicely {
    my $msg = shift;
    
    print $msg;
    exit;
}
