#!/usr/bin/perl

use strict;

my $currentcluster;
my $currentlen;
my $currentname;
my $outputfasta = $ARGV[0];
my $index = $ARGV[1];
my $coverage = $ARGV[2];

my $outputfh;
open $outputfh, "+>$outputfasta" or die "Can't open output file $outputfasta\n";

my $idlookup = {};
my $fastafile;
open $fastafile, "$index" or die "Can't open fasta file";
while(my $line=<$fastafile>){
    if($line =~ /^>/){
	chomp $line;
	my($id) = ($line =~ /^>(\S+)/);
	$idlookup->{$id}=1;
    }
}


while(my $line=<STDIN>){
    chomp $line;
    my($clustername,$seq,$length) = split(/\t/,$line);
    if($clustername ne $currentcluster){
	&printrep($seq,$index,$outputfh);
	$currentcluster=$clustername;
	$currentlen=$length;
	$currentname=$seq;
	$idlookup->{$seq}=0;
	$idlookup->{$currentname}=0;
	print "$currentname\t$seq\n";
    }
    else{
	if($length/$currentlen >= $coverage){
	    $idlookup->{$seq}=0;
	    print "$currentname\t$seq\n";
	}
	else{
	    &printrep($seq,$index,$outputfh);
	    $currentcluster=$clustername;
	    $currentlen=$length;
	    $currentname=$seq;
	    $idlookup->{$seq}=0;
	    $idlookup->{$currentname}=0;
	    print "$currentname\t$seq\n";
	}
    }
}

&printremainders($outputfh,$idlookup);

sub printrep{
    my($seq,$index,$outputfh) = @_;
    #my($db,$acc) = ($seq =~ /(\w+)\|(\w+)/);
    #print seq in fasta using index
    #print $outputfh `xdget -T$db -p $index $acc`;
    print $outputfh ">$seq\n";

}

sub printremainders{
    my($outputfh,$idlookup) = @_;
    foreach my $id (keys %$idlookup){
	print $outputfh ">$id\n" if($idlookup->{$id}==1);
    }
}
