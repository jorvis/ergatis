#!/usr/local/bin/perl

use lib("../");
use strict;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);


my %options = ();
my $results = GetOptions (\%options,  'pegene|g=s', 'pematch|m=s', 'verbose|v', 'pebin|b=s', 'help|h' );


###-------------PROCESSING COMMAND LINE OPTIONS-------------###

my $pegene     = $options{'pegene'};
my $pematch    = $options{'pematch'};
my $verbose    = $options{'verbose'};
my $pebin      = $options{'pebin'} || "/usr/local/devel/ANNOTATION/shared/bin/linux/peffect";
my $QUERYPRINT;
my $DEBUG = $options{'DEBUG'} || 0;

if(!$pegene or !$pematch or exists($options{'help'})) {
    &print_usage();
}

###-------------------------------------------------------###

if(! -s $pegene) {
    print STDERR "Unable to find $pegene.  Aborting.Now..\n";
    exit 2;
}
if(! -s $pematch) {
    print STDERR "Unable to find $pematch.  Aborting...\n";
    exit 3;
}



#system("$pebin  -w 10 -g -50 -r -100 -m 4 -o 3 -f $pegene < $pematch") and die "Unable to execute\n";#
my $output = qx($pebin  -w 10 -g -50 -r -100 -m 4 -o 3 -f $pegene < $pematch 2>&1);
if($?) {
    print STDERR "Unable to execute $pebin\n";
    exit 5;
}else {
    print "$output";
}



sub print_usage {


    print STDERR "SAMPLE USAGE:  run_pe.pl -g 1vs5.xml -m 1vs5.match.xml > pe.out\n";
    print STDERR "  --pegene(-g)  = output of db2xml.pl \n";
    print STDERR "  --pematch(-m) = output of dbmatch2xml.pl\n";
    print STDERR "  --pebin       = full path to pebin (optional)\n";
    print STDERR "  --help = This help message.\n";
    exit 1;

}









