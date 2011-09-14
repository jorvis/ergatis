#!/usr/bin/env perl

=head1 NAME

genbank2fasta.pl - Takes genbank files and create fasta

=head1 SYNOPSIS

 USAGE: genbank2fasta.pl
       --input_file=/path/to/some/input.file
       --output=/path/to/transterm.file
     [ --log=/path/to/file.log
       --debug=3
       --help
     ]

=head1 OPTIONS

B<--input_file,-i>
	Should be an input fasta file

B<--output_file,-o>
	Will print out the molecule fasta.

B<--help,-h>
    Print this message

=head1  DESCRIPTION

	If no sequence information is in the genbank file the script will fail
 
=head1  INPUT
	Genbank record

=head1 OUTPUT
	Fasta output.
	Header:
	
	TODO

=head1  CONTACT

    Kevin Galens
    kgalens@gmail.com

=cut

use strict;
use warnings;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Pod::Usage;
use Bio::SeqIO;

my %options;
my $results = GetOptions (\%options,
                         "input_file|i=s",
                         "output_file|o=s",
                         "help|h"
                          );

&check_options(\%options);


my $in = Bio::SeqIO->new( -file => $options{'input_file'},
							-format => 'genbank' );
my $out = Bio::SeqIO->new( -file => ">$options{'output_file'}",
						   -format => 'fasta' );

while( my $seq = $in->next_seq() ) {
	$out->write_seq($seq) ;
}

sub check_options {

   my $opts = shift;

   if( $opts->{'help'} ) {
       &_pod;
   }

   foreach my $req ( qw(input_file output_file) ) {
       die("Option $req is required") unless( $opts->{$req} );
   }
}

sub _pod {
    pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDERR} );
}
