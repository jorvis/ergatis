#!/usr/local/bin/perl
#-----------------------------------------------------------------------
# util:     watchwf.pl
# author:   Jay Sundaram
# date:     2004-03-27
#
# purpose:  Alternative lightweight method for monitoring a Workflow
# 
#-----------------------------------------------------------------------
use Data::Dumper;
use strict;

my $string = $ARGV[0];

if (!defined($string)){
    print STDERR "You must specify the Workflow instance file\n";
    exit(1);
}

if (!-e $string){
    print STDERR "Workflow instance file '$string' does not exist\n";
    exit(1);
}

if (!-T $string){
    print STDERR "Workflow instance file '$string' is not a text file\n";
    exit(1);
}


my $sleep  = $ARGV[1];
if (!defined($sleep)){
    $sleep = 10;
    print STDERR "sleep set to '$sleep'\n";
}

print STDERR "Remember to kill this job when you're done (process id $$)\n";

my $checks =  {
	       states   => "grep state $string | wc -l",
	       complete => "grep state $string | grep \">complete\" | wc -l",
	       running  => "grep state $string | grep \">running\" | wc -l",
	       failed   => "grep state $string | grep \">failed\" | wc -l",
	       pending  => "grep state $string | grep \">pending\" | wc -l",
	       incomplete =>  "grep state $string | grep \">incomplete\" | wc -l"
	   }; 

while (1){
    sleep $sleep;
    
    my $date = qx{date};
    chomp $date;
    print "$date Monitoring $string (with process id $$)\n";
    foreach my $key (sort keys %$checks){
	my $grep = $checks->{$key};
	
	my $stat = qx{$grep};
	my $num;
	if ($stat =~ /(\d+)/){
	    $num =$1;
	}
	else{
	    die "wc failed...";
	}
	printf "%-11s    %-2d\n", $key, $num;
    }
    print "\n";
}
