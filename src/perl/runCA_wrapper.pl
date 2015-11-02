#!/usr/bin/perl

## Will run the runCA executable with the passed in options.
## This will take in a list of files and call runCA using all 
## of the files as input

use strict;
use warnings;
use Getopt::Long;

my %options;
my $results = GetOptions (\%options, 
                          'input_list|i=s',
                          'runca_opts|r=s',
                          'runca_bin|b=s' );

foreach my $req ( qw(input_list runca_opts runca_bin) ) {
    die("Option $req is required") unless( $options{$req} );
}

open(IN, "< $options{'input_list'}") or die("Couldn't open $options{'input_list'} $!");
chomp( my @files = <IN> );
close(IN);
my $input_string = join(" ", @files);

my $cmd = $options{'runca_bin'};
$cmd .= " ".$options{'runca_opts'};
$cmd .= " ".$input_string;

open( OUT, "$cmd |") or die("Could not open $cmd");
map { print $_; } <OUT>;
close(OUT);


