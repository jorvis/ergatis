#!/usr/bin/perl
#Input from mummer
#Output 
#qry\tref\trefstart\trefend\tqrystart\tqryend
#Note: querymatchlen == refmatchlen 

use strict;

my $qryid;
my $lengths={};
while(my $line=<STDIN>){
    chomp $line;
    if($line =~ /^\>\s*(\S+)\s+Len\s+=\s+(\d+)/){
	$qryid = $1;
	$lengths->{$1} = $2;
    }
    elsif($line =~ />\s*(\S+)/){
	$qryid = $1;
    }
    elsif($line !~ /^\#/){
	$line =~ s/^\s+//g;
	my($refid,$startref,$startqry,$len) = split(/\s+/,$line);
	print "$refid\t$qryid\t$startref\t",$startref+$len-1,"\t$startqry\t",$startqry+$len-1,"\t",$lengths->{$qryid},"\n";
    }
}
