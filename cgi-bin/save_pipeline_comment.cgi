#!/usr/local/bin/perl -w

use strict;
use CGI;
use CGI::Carp qw(fatalsToBrowser);

my $q = new CGI;
my $comment = $q->param('pipeline_comment') || carp("pipeline_comment is a required argument");
my $file = $q->param('pipeline_comment_file') || carp("pipeline_comment_file is a required argument");

print $q->header( -type => 'text/plain' );

open( my $ofh, ">$file" ) || carp("failed to create comment file: $!");

## if the entire comment is whitespace, clear it
if ( $comment =~ /^\s*$/ ) {
    $comment = '';
}

print $ofh $comment;

## the ajax call might do something with it too
print "$comment";
