#!/usr/bin/perl

use strict;

my $file = shift @ARGV;
open FILE,"$file" or die "Can't open $file";
my @files = <FILE>;
chomp @files;

my $cmd = "sort -m ".join(' ',@ARGV)." ".join(' ',@files);
print STDERR "Running $cmd\n";
print $cmd;
