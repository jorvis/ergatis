#!/usr/bin/perl -w
#mergelines.pl

use strict;

my $count = 0;
while(my $line = <STDIN>){
    chomp($line);
    print $line;
    $count = ($count + 1) % 4;
    
    if($count == 0){
        print "\n"; 
    }else{
        print "\t";
    }
}

