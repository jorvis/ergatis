#!/usr/bin/perl

use strict;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use MLDBM 'DB_File';

my %info;

my $db = shift || die "pass me a database\n";

## create the tied hash
tie(%info, 'MLDBM', $db ) || die "failed to tie: $!";

for my $acc ( keys %info ) {
    print "$acc\n";
    
    for my $key ( keys %{$info{$acc}} ) {
        print "\t$key\t=>\t$info{$acc}{$key}\n";
    }
    
    #print "\n\nenter to see another ... ";
    #<STDIN>;
    print "\n\n";
}

untie(%info);

exit(0);

