#!/usr/bin/perl -w

use strict;
use CGI;
use CGI::Carp qw(fatalsToBrowser);

my $q = new CGI;
my $comment = $q->param('comment') || carp("pipeline_comment is a required argument");
my $file = $q->param('comment_file') || carp("pipeline_comment_file is a required argument");

print $q->header( -type => 'text/plain' );

## make sure the directory exists
$file =~ m|^(.+)/|;
my $dir = $1;

if ( ! -e $dir ) {
    mkdir($dir) || die "can't create directory for pipeline comment: $!";
}

open( my $ofh, ">$file" ) || carp("failed to create comment file: $!");

## if the entire comment is whitespace, clear it
if ( $comment =~ /^\s*$/ ) {
    $comment = '';
}

print $ofh $comment;

## the ajax call might do something with it too
print "$comment";
