#!/usr/bin/perl -w

use strict;
use CGI;
use CGI::Carp qw(fatalsToBrowser);

my $q = new CGI;
my $path = $q->param('path') || '';

my $layout_path = $path . '/pipeline.layout';

print $q->header( -type => 'text/plain' );

if ( -e $layout_path ) {
    open(my $ifh, $layout_path) || croak( "failed to open layout file $layout_path: $!" );
    
    print <$ifh>;

} else {
    croak( "$layout_path not found" );
}
