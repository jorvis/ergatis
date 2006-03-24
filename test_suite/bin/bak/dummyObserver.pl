#!/usr/bin/perl

#use warnings;
#use strict;
use Cwd qw(realpath);
use File::Basename;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);

my %param=();
GetOptions(\%param, "ID=s", "name=s", "file=s", "host=s", "event=s", "message=s","time=s","props=s","retval=s");

my $exec_path = dirname(realpath($0));
open (LOG, ">>$exec_path/dummyObserver.log") || die "Failed opening log for writing\n";
#print LOG join(" ", @ARGV)."\n";
print LOG $param{'ID'}."\t".$param{'file'}."\t".$param{'event'}."\t".$param{'time'}."\n";
close LOG;
exit(0);
