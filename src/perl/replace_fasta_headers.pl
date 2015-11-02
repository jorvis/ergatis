#!/usr/bin/perl

=head1  NAME 

replace_fasta_headers.pl - Replace fasta headers with a short code (prot_#) and create a mapping file

=head1 SYNOPSIS

USAGE: replace_fasta_headers.pl --input=/path/to/input --map_file=/path/to/map_file --replace --output=/path/to/output

=head1 OPTIONS

B<--input,-i> 
    Input file file 

B<--debug,-d> 
    Debug level.  Use a large number to turn on verbose debugging. 

B<--log,-l> 
    Log file

B<--map_file,-m> 
    path to the map file, whether it has been created or not. 

B<--replace,-r> 
    whether or not to replace the keys in the input with the values in the map file. Not having this option
    will generate the map file (assuming your input is a fasta file).

B<--output,-o> 
    path to the output file. This is used to specify where the output fasta will be written. 

B<--help,-h> 
    This help message

=head1   DESCRIPTION

This script will take a fasta file and replace the headers with checksums. It will also generate a mapping file
that can be used to replace the checksums later on. 

=head1 INPUT

=head1 OUTPUT


=head1 CONTACT

David Riley
driley@som.umaryland.edu

=cut

use strict;
use Log::Log4perl qw(get_logger);
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use Pod::Usage;

my %options = ();
my $results = GetOptions (\%options, 
              'input|i=s',
              'map_file|m=s',
              'output|o=s',
              'debug|d=s',
              'replace=s',
              'log|l=s',
              'do_nothing=s',
              'help|h') || pod2usage();

# display documentation
if( $options{'help'} ){
    pod2usage( {-exitval=>0, -verbose => 2, -output => \*STDOUT} );
}

exit if($options{'do_nothing'});

if($options{'replace'} ) {
    &replace_checksums();
}
else {
    &generate_checksums();
}

sub generate_checksums {
    
    open IN, "<$options{'input'}" or die "Unable to open input file $options{'input'}";
    open MAP, ">$options{'map_file'}" or die "Unable to open map_file $options{'map_file'}";
    open OUT, ">$options{'output'}" or die "Unable to open output fasta file $options{'output'}";
    my $count = 1;
    while(<IN>) {
        if(/^>\s*(\S+)/) {
            print MAP "prot_$count\t$1\n";
            print OUT ">prot_$count\n";
            $count++;
        }
        else {
            print OUT $_;
        }
    }
    close IN;
    close OUT;
    close MAP;
}

sub replace_checksums {
    open IN, "<$options{'map_file'}" or die "Unable to open map file $options{map_file}\n";

    while(<IN>) {
        chomp;
        my ($key, $value) = split(/\t/, $_);
        print "replacing $key with $value in $options{'input'}\n";
        `perl -pi -e 's[$key][$value]g' $options{'input'}`;
    }
    close IN;
}
