#!/usr/bin/env perl

=head1 NAME

create_filter_fastq_file_iterator_list.pl

=head1 SYNOPSIS

 USAGE: create_filter_fastq_file_iterator_list.pl
       --input_map=/path/to/input_map.txt
       --output_iter_list=/path/to/i1.list

=head1 OPTIONS

B<--input_map>
    Input file should be tab delimited.
    contigs.fa    fastq1[,fastq2]    (Sequences|Read List)

B<--output_iter_list>

=head1  CONTACT

    Kevin Galens
    kgalens@gmail.com

=cut

use strict;
use warnings;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use File::Basename;

my %options;
my $results = GetOptions (\%options,"input_map=s","output_iter_list=s");
&check_options(\%options);

open(IN, "< $options{'input_map'}") or die("Can't open $options{'input_map'}: $!");
open(OUT, "> $options{'output_iter_list'}") or die("Can't open $options{'output_iter_list'} for writing: $!");

print OUT "\$;BASENAME\$;\t\$;READS\$;\t\$;READ_LIST\$;\n";

while( my $line = <IN> ) {
    chomp($line);
    my ($fa, $reads, $read_list) = split(/\t/, $line);
    my $unique_id = &get_id( $reads );
    print OUT "$unique_id\t$reads\t$read_list\n";
}

close(IN);
close(OUT);

sub get_id {
    my ($reads) = @_;
    my @files = split(/,/,$reads);
    basename( $files[0], qw(.fq .fastq .fq.gz .fastq.gz) );
}

sub check_options {
   my $opts = shift;

   foreach my $req ( qw(input_map output_iter_list) ) {
       die("Option $req is required") unless( $opts->{$req} );
   }
}
