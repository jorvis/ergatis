#!/usr/local/bin/perl

use strict;

use CGI qw(:standard);

use File::Basename;

my $file = param('file');
my $removedir = param('removedir');

if($removedir){
    my $dirname = dirname($file);
    
    if($dirname =~ /\d+$/){
	`rm -rf $dirname`;
	print header();
	print "<html><body>Removed directory $dirname<br><a href='javascript:window.parent.location.reload();window.close()'>[close]</a></body></html>";
	exit;
    }
    else{
	`rm -rf $file`;
	print header();
	print "<html><body>Removed file $file<br><a href='javascript:window.parent.location.reload();window.close()'>[close]</a></body></html>";
	exit;
    }
}
else{
    `rm -rf $file`;
    print header();
    print "<html><body>Removed file $file<br><a href='javascript:window.parent.location.reload();window.close()'>[close]</a></body></html>";
    exit;
}


