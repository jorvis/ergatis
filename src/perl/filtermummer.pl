#!/usr/bin/perl

use strict;

use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);

my %options = ();
my $results = GetOptions (\%options, 
                          'reference|r=s',
                          'query|q=s');

while(my $file=<STDIN>){
    chomp $file;
    if(-e $file){
	open FILE,"$file" or die "Can't open $file";
	while(my $line=<FILE>){
	    chomp $line;
	    my @x=split(/\s+/,$line);
	    if(exists $options{'query'}){
		my $qrycov= ($x[5]-$x[4]+1)/$x[6];
		print "$line $qrycov\n" if($x[0] ne $x[1] && $qrycov >= $options{'query'} && $x[5] eq $x[6]);
	    }
	    if(exists $options{'reference'}){
		my $refcov = ($x[3]-$x[2]+1)/$x[8];
		print "$line $refcov\n" if($refcov >= $options{'reference'} && $x[3] eq $x[8]);
	    }
	}
    }
}
