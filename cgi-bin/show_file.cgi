#!/usr/local/bin/perl

use strict;

use CGI qw(:standard);

my $file = param('file');

print header();

print "<html><body><pre>";

open FILE, $file or print "Can't open file $file\n";
while (my $line=<FILE>){
    $line = &link_ip($line);
    print $line;
}
close FILE;

print "</pre></body></html>";

sub link_ip{
    my($line) = @_;
    if($line =~ /\d{1,3}\.\d{1,3}\.\d{1,3}.\d{1,3}/){
	my(@ips) = ($line =~ /(\d{1,3}\.\d{1,3}\.\d{1,3}.\d{1,3})/g);
	foreach my $ip (@ips){
	    print STDERR "Looking for $ip\n";
	    my $host = `host $ip`;
	    if($host !~ /not found/){
		chomp $host;
		$host =~ s/.*pointer //;
		my $url = "http:\/\/intranet\/sysadmin\/bb\/html\/$host"."cpu.html";
		$line =~ s/$ip/<a href="$url">$host<\/a>/;
	    }
	}
    }
    return $line;
}
