#! /local/perl/bin/perl

my $str =  "cat $ARGV[0]".' | '."$ARGV[1]".'>'."$ARGV[2]";
system $str;
