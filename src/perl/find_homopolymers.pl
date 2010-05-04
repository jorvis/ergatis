#!/usr/bin/perl

eval 'exec /usr/bin/perl  -S $0 ${1+"$@"}'
    if 0; # not running under some shell

=head1  NAME 

find_homopolymers.pl - Find homopolymeric tracts in a sequence file.

=head1 SYNOPSIS

USAGE: find_homopolymers.pl --input=/path/to/somefile.fsa --at_min_length=8 --gc_min_lenght=5 [--output=/path/to/somefile.out] [--help]
      
=head1 OPTIONS

=over 8

=item B<--input,-i> 

Input fasta file.

=item B<--output,-o> 

Optional. Output file (if not provided the results will print to STDOUT)

=item B<--at_min_length,-at> 

Minimum length of an A/T repeat

=item B<--gc_min_length,-gc> 

Minimum length of an G/C repeat

=item B<--help,-h> 

This help message

=back

=head1   DESCRIPTION

This script finds homopolymeric tracts in a DNA sequence

=head1 INPUT

Input sequence in fasta format.

=head1 OUTPUT

List of coordinates of homopolymers and their sequence.

=head1 CONTACT

    David Riley
    driley@som.umaryland.edu

=cut

use strict;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use Pod::Usage;


my %options = ();
my $results = GetOptions (\%options, 
              'input|i=s',
              'output|o=s',
              'at_min_length|at=s',
              'gc_min_length|gc=s',
              'help|h');

# display documentation
if( $options{'help'} ){
    pod2usage( {-exitval=>0, -verbose => 2, -output => \*STDOUT} );
}

## make sure all passed options are peachy
&check_parameters(\%options);

my $at_min_length=$options{'at_min_length'};
my $gc_min_length=$options{'gc_min_length'};
my $seqfile = $options{'input'};
my $output;
if($options{'output'}) {
open $output, ">$options{'output'}" or die "Unable to open $options{'output'}\n";
}
else {$output = *STDOUT};

open IN, "<$seqfile" or die "Unable to open input file $seqfile\n";;
my @seq;
my $id = '';
while(<IN>) {
    chomp;
    if($_!~/^>/) {
        push(@seq, $_);
    }
    else {
        $id = $_;
    }
}

my $seqstring = join('',@seq);
print $output "$id (".length($seqstring)."):\n";
while ($seqstring =~ /([Aa]{$at_min_length,}|[Tt]{$at_min_length,}|Gg]{$gc_min_length,}|[Cc]{$gc_min_length,})/g) {
    print $output length($`)."\t".(length($`)+length($&))."\t".$1."\n";
}

sub check_parameters {
    
    ## make sure input file exists
    if (! -e $options{'input'}) { print STDERR "Input file invalid\n"; pod2usage( {-exitval=>0, -verbose => 0, -output => \*STDOUT})};

    if(!$options{'at_min_length'}) { pod2usage( {-exitval=>0, -verbose => 0, -output => \*STDOUT} )};
    if(!$options{'gc_min_length'}) { pod2usage( {-exitval=>0, -verbose => 0, -output => \*STDOUT} )};

    return 1;
}
