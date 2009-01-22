#!/usr/bin/perl

use strict;

my $fastanames = {};
my $fastalen = {};
my $names=0;
my $len=0;
print "strict graph mygraph{\n";

while(my $line = <STDIN>){
    chomp $line;
    my @x=split(/\s+/,$line);
    if(! exists $fastanames->{$x[0]}){
	$fastanames->{$x[0]}=++$names;
	print "$names [fastaname=\"$x[0]\",fastalength=\"$x[8]\"];\n";
    }
    if( ! exists $fastanames->{$x[1]}){
	$fastanames->{$x[1]} =++$names; 
	print "$names [fastaname=\"$x[1]\",fastalength=\"$x[6]\"];\n";
    } 
    my $refcov = $x[7];
    my $qrycov = $x[9];
    if($x[7]<1){
	($refcov) = ($x[7] =~ /(.\d)/);
    }
    if($x[9]<1){
	($qrycov) = ($x[9] =~ /(.\d)/);
    }
    print "$fastanames->{$x[0]} -- $fastanames->{$x[1]} [label=\"$refcov,$qrycov\"]\n";
}
print "}\n";
