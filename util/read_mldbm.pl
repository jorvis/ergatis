#!/usr/bin/perl

use strict;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use MLDBM 'DB_File';

my %info;

my $db = shift || die "pass me a database\n";
my $accession = shift || '';

## create the tied hash
tie(%info, 'MLDBM', $db ) || die "failed to tie: $!";

if ($accession) {
    my $entry = $info{$accession} || die "failed to find accession: $accession";
    print "$accession\n";
    for my $key ( keys %$entry ) {
        print "\t$key\t=>\t$$entry{$key}\n";
    }
    
} else {
    for my $acc ( keys %info ) {
        print "$acc\n";

        for my $key ( keys %{$info{$acc}} ) {
            print "\t$key\t=>\t$info{$acc}{$key}\n";
        }
    }
}
print "\n\n";


untie(%info);

exit(0);

