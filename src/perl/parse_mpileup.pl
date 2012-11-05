#!/usr/bin/env perl

=head1 NAME

parse_mpileup.pl - 

=head1 SYNOPSIS

 USAGE: parse_mpileup.pl
       --merged_table=input
       --output=/path/to/transterm.file
     [ --log=/path/to/file.log
       --debug=3
       --help
     ]

=head1 OPTIONS

=head1  DESCRIPTION

 DESCRIPTION
 
=head1  INPUT
    Describe the input

=head1 OUTPUT
    Describe the output

=head1  CONTACT

    Kevin Galens
    kgalens@gmail.com

=cut

use strict;
use warnings;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use List::Util qw(sum);
use Data::Dumper;
use Pod::Usage;

############# GLOBALS AND CONSTANTS ################
####################################################

my %options;
my $results = GetOptions (\%options,
			  "mpileup|p=s",
			  "mapped_positions|m=s",
			  "query|q=s",
			  "version|v=s",
			  "output|o=s",
                          );

&check_options(\%options);

print "Parsing mapped positions\n";
my $mapped_positions = &parse_query_positions( $options{'mapped_positions'}, $options{'query'} );

print "Parsing mpileup\n";
my $qual_cov = &parse_mpileup( $options{'mpileup'}, $mapped_positions, $options{'version'} );

print "Printing output\n";
&print_output( $qual_cov, $options{'output'} );

sub print_output {
  my( $data, $file ) = @_;

  open(OUT, "> $file") or die("Can't open $file for writing: $!");
  foreach my $refmol ( sort { $a cmp $b } keys %{$data} ) {
	foreach my $refpos( sort { $a <=> $b } keys %{$data->{$refmol}} ) {
	  my %data = ( 'qmol' => [], 'qpos' => [], 'cov' => [], 'qual' => [] );
	  foreach my $entry ( @{$data->{$refmol}->{$refpos}} ) {
		push(@{$data{'qmol'}}, $entry->[0] );
		push(@{$data{'qpos'}}, $entry->[1] );
		push(@{$data{'cov'}}, $entry->[2] );
		push(@{$data{'qual'}}, $entry->[3] );
	  }
	  print OUT "$refmol\t$refpos";
	  print OUT "\t".join("/", @{$data{$_}}) for( qw(qmol qpos cov qual) );
	  print OUT "\n";
	}
  }
  close(OUT);
  
}

sub parse_mpileup { 
  my ($mpileup, $query_positions, $version) = @_;

  my $results = {};

  open(IN, "<$mpileup") or die("Can't open $mpileup: $!");
  while(my $line = <IN>) {
	chomp( $line );
	my @c = split(/\t/, $line);

	# skip any row that we don't care about [those rows where we don't have a SNP]
	next unless( exists( $query_positions->{$c[0]}->{$c[1]} ) );

	# coverage is easy
	my $coverage = $c[3];

	# and quality
	my $av_qual = &get_average_quality($c[5], $version);

	foreach my $arr ( @{$query_positions->{$c[0]}->{$c[1]}} ) {
	    my ($refmol, $refpos, $query) = @{$arr};
	    $results->{$refmol} = {} unless( exists( $results->{$refmol} ) );
	    $results->{$refmol}->{$refpos} = [] unless( exists( $results->{$refmol}->{$refpos} ) );
	    push(@{$results->{$refmol}->{$refpos}}, [$c[0], $c[1], $coverage, $av_qual]);
	}

 }
  close(IN);

  return $results;
}

sub get_average_quality {
  my ($qual_string, $version) = @_;
  my $offset = 0;
  if( $version eq 'new' ) {
	$offset = 33;
  } elsif( $version eq 'old' ) {
	$offset = 64;
  } else {
	die("Could not understand version $version.");
  }

  ## Different for perl version 5.8 than for version 5.10
  my @unpacked;
  if( $] < 5.009 ) {
	@unpacked =  unpack("C*", $qual_string);
  } else {
	@unpacked =  unpack("W*", $qual_string);
  }

  my $av_quality = sprintf( "%.2f", sum( map { $_ -= $offset } @unpacked )/scalar(@unpacked) );

  return $av_quality;
}

sub parse_query_positions {
  my ($file, $query ) = @_;

  my $positions_by_query = {};

  open( IN, "< $file") or die("Can't open $file: $!");
  chomp( my $h = <IN> );
  my @header = split( /\t/, $h );

  # Find the index of the query we are dealing with
  my @index = grep { $header[$_] eq $query } 0..$#header;
  die("Could not find query $query in header of mapped positions file. [@index]") unless( @index == 1 );
  my $i = $index[0];  

  while( my $line = <IN> ) {
	chomp( $line );
	my @c = split(/\t/, $line);

	# Get the contig id and query position
	my ($contig, $pos) = ($c[$i-1], $c[$i]);
	die("Could not parse contig from line [$line]") unless( $contig );
	die("Could not parse pos from line [$line]") unless( $pos );

	# Sometimes we have multiple mappings.
	my @contigs = split(m|/|, $contig);
	my @positions = split(m|/|, $pos);

	# Make sure we had same number of contig names as positions
	unless( @contigs == @positions ) {
	  print Dumper( \@contigs );
	  print Dumper( \@positions );
	  die("Didn't find the same number of contigs as positions");
	}

	# cycle through each contig id and position and store the mapping
	for ( my $j = 0; $j < @contigs; $j++ ) {
		
	  # I assume that multiple SNPs don't map to the same contig location. Just
	  # in case that's not true, I'm checking for it here.
	  if ( !exists( $positions_by_query->{$contigs[$j]}->{$positions[$j]} ) ) {
	      $positions_by_query->{$contigs[$j]}->{$positions[$j]} = [];
	  }
	  
	  push(@{$positions_by_query->{$contigs[$j]}->{$positions[$j]}}, [$c[0], $c[1], $query]);
	}

  }
  close(IN);

  $positions_by_query;
}


sub check_options {
   my $opts = shift;
   
   foreach my $req ( qw(version query mpileup mapped_positions output) ) {
	 die("Option $req is required") unless( $opts->{$req} );
   }
}
