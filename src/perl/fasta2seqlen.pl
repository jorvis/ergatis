#!/usr/bin/perl
#Input FASTA
#Output 
#FASTAheaderid\tseqlen\tuid



use strict;

my $headerid;
my $len=0;
my $uid=0;
while(my $line=<STDIN>){
    chomp $line;
    if($line =~ /^\>\s*(\S+)/){
	print "$headerid\t$len\n" if(defined $headerid && $len>0);
	$headerid = $1;
	$len=0;
	$uid++;
    }
    elsif($line =~ /^\S+/){
	$line =~ s/\s+//g;
	$len += length($line);
    }
}
print "$headerid\t$len\n" if(defined $headerid && $len>0);
