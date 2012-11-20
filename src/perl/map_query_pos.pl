#!/usr/bin/env perl
use strict;
use warnings;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use File::Basename;
use SNP::MergedTable;

use Data::Dumper;

my %options = &check_options();

## Holds a map of query accessions to fasta files
my %blasted_queries; 

## %query_map contains a map of contig ids to query accession
my %query_map = &generate_query_map( $options{'query_acc_map'} );

## Holds the query names from the merged table
my @queries;

## ref_positions: $ref_positions{$ref_molecule}->{$ref_pos} = [$refbase,$properties]
my %ref_positions = &parse_merged_table( $options{'merged_table'} );

my %query_positions;

foreach my $refmol ( keys %ref_positions ) {
  print "$refmol: ".scalar(keys %{$ref_positions{$refmol}})."\n";
}

&parse_blast_list( $options{'raw_blast_list'} );

open(my $ofh, "> $options{'output_file'}") or die("Can't open $options{'output_file'} for writing: $!");
&write_to_outfile( $ofh );
close($ofh);

sub check_options {
    my %options;
    my $results = GetOptions (\%options,
			      'query_acc_map|q=s',
			      'raw_blast_list|r=s',
			      'merged_table|m=s',
			      'output_file|o=s'
                              );
    
    foreach my $req ( qw(query_acc_map raw_blast_list merged_table output_file) ) {
        die("Option $req is required") unless( $options{$req} );
    }
    return %options;
}

sub write_to_outfile {
  my ($ofh) = @_;

  my @bqueries = sort keys %blasted_queries;

  my @header = ("refmol", "refpos", map { ("", $_) } @bqueries);
  print $ofh join("\t", @header)."\n";
  
  foreach my $refmol ( sort { $a cmp $b } keys %ref_positions ) {
	foreach my $refpos( sort { $a <=> $b } keys %{$ref_positions{$refmol}} ) {
	  my @c;
	  push(@c, ($refmol, $refpos) );

	  foreach my $query ( @bqueries ) {
	      if( exists( $query_positions{$refmol}->{$refpos}->{$query} ) ) {

		  my $maps = $query_positions{$refmol}->{$refpos}->{$query};
		  my $best_eval = 1;
		  my @best_mappings;
		  foreach my $map ( @{$maps} ) {
		      if( $map->[2] < $best_eval ) {
			  @best_mappings = ($map);
			  $best_eval = $map->[2];
		      } elsif( $map->[2] == $best_eval ) {
			  push(@best_mappings, $map);
		      }
		  }

		  push(@c, (join("/",map{ $_->[0] } @best_mappings), join("/",map{ $_->[1] } @best_mappings)) );
		  
	      } else {
		  print "No mapping for $refmol:$refpos on query $query\n";
		  #die("Could not find mapping for $query [$refmol, $refpos]");
	      }
	  }
	  print $ofh join("\t", @c)."\n";
      }
    }

}

sub parse_merged_table {
  my ($merged_table) = @_;
  my %retval;
  print "Parsing merged table:\n";
  my $mtable = SNP::MergedTable::parse( $merged_table );

  @queries = $mtable->queries;

  my @rows = $mtable->get_rows();

  my $total = scalar(@rows);
  my $count = 0;
  print "$count / $total\r";
  foreach my $row ( @rows ) {
	$retval{$row->molecule} = {} unless( exists( $retval{$row->molecule} ) );
	$retval{$row->molecule}->{$row->refpos} = [$row->refbase,$row->properties];
	$count++;
	print "$count / $total\r";
  }
  print "\nDone.\n";
  undef $mtable;
  return %retval;
}

sub parse_blast_list {
  my ($list) = @_;

  print "Parsing BLAST list:\n";
  open(IN, "< $list") or die("Can't open $list: $!");
  chomp( my @files = <IN> );
  close(IN);

  my $total = scalar(@files);
  my $count = 0;

  my $snps = 0;
  foreach my $blast ( @files ) {
      print "$blast [ $count / $total ]\n";
      $snps += &parse_blast_file( $blast );
      $count++;
      #print "$count / $total\r";
  }
  
  print "\nDone. Found $snps SNPS\n";
}

sub parse_blast_file {
  my ($file) = @_;

  open(my $in, "< $file") or die("Can't open $file: $!");

  my $count = 0;

  my $flag = 0;
  my ($refmol, $pos);
  while( my $line = <$in> ) {
	chomp( $line );
	if( $flag && $line =~ /^>(\S+)/ ) {
	  my $id = $1;
	  
	  # does the id exist? can we map it?
	  unless( exists( $query_map{ $id } ) ) {
	      die("Couldn't map blast subject: $id");
	  }

	  my ($subject_pos,$eval) = &parse_hit( $in, $id );
	  &store_position( $id, $subject_pos, $eval, $refmol, $pos );
	  $count++;
	}

	if( $line =~ /^Query=\s*(\S+)_SNP_(\d+)/ ) {
	  ($refmol, $pos) = ($1,$2);
	  $flag = (exists( $ref_positions{$refmol}->{$pos} )) ? 1 : 0;
	}
  }
  close($in);
  
  return $count;
}

sub store_position {
  my ($id, $pos, $eval, $refmol, $refpos) = @_;
  
  # Get the query ID from the map
  my $acc = $query_map{ $id };
  die("Could not find contig with id: $id") unless( $acc );

  $query_positions{$refmol} = {} unless( exists( $query_positions{$refmol} ) );
  $query_positions{$refmol}->{$refpos} = {} unless( exists( $query_positions{$refmol}->{$refpos} ) );
  $query_positions{$refmol}->{$refpos}->{$acc} = [] unless( $query_positions{$refmol}->{$refpos}->{$acc} );
  push(@{$query_positions{$refmol}->{$refpos}->{$acc}}, [$id, $pos, $eval]);
}

sub parse_hit {
  my ($fh) = @_;

  my $qstart;
  my $sstart;
  my $send;
  my $expect;
 
  my $count = 0;
  while( my $line = <$fh> ) {
	chomp( $line );
	if( $line =~ /^Query:\s*(\d+)/ ) {
	  $qstart = $1;
	} elsif( $line =~ /^Sbjct:\s*(\d+)\s*\S+\s*(\d+)/ ) {
	  $sstart = $1;
	  $send = $2; 
	  die("Could not parse subject start and end [$sstart - $send]") 
		unless( defined( $sstart ) && defined( $send ) );
	  last;
	} elsif( $line =~ /Expect\s*=\s*(\S+)/ ) {
	  $expect = $1;
	}
	$count++;
	last if( $count > 20 );
  }

  die("Could not parse qstart or sstart or subject_end") unless( $qstart && $sstart && $send );
  die("Could not parse expect value") unless( $expect );

  my $reverse = ( $send < $sstart ) ? 1 : 0;
  my $offset = 20 - ($qstart - 1);

  my $pos = $sstart + $offset;
  $pos = $sstart - $offset if( $reverse );

  return ($pos, $expect);
}


sub generate_query_map {
  my ($list) = @_;
  my %retval;

  print "Generating query header map:\n";
  open(IN, "< $list") or die("Can't open $list: $!");
  chomp( my @lines = <IN> );
  close(IN);

  my $acc_fsa_map = {};
  foreach my $l ( @lines ) {
      my @t = split(/\t/, $l); 
      die("Didn't parse 2 columns from line [$l]. Parsed ".scalar(@t)) unless( @t == 2 );
      $acc_fsa_map->{$t[0]} = $t[1];
  }

  my $total = scalar(keys %{$acc_fsa_map});
  my $count = 0;
  print "$count / $total\r";

  foreach my $query_acc ( keys %{$acc_fsa_map} ) {
      my $f = $acc_fsa_map->{$query_acc};
      $blasted_queries{$query_acc} = $f;

      open(IN, "< $f" ) or die("Can't open $f: $!");
      my @headers = grep { /^>/ } <IN>;
      close(IN);

      foreach my $h ( @headers ) {
	  my $key = $1 if( $h =~ /^>(\S+)/ );
	  die("Could not parse key from header $h") unless( $key );
	  $key = $1 if( $key =~ /^gi\|\d+\|([^\|]+\|[^\|]+\|)/ );
	  $retval{$key} = $query_acc;
      }

      $count++;
      print "$count / $total\r";
  }
  print "\nDone.\n";

  return %retval;
}


=head1 NAME

map_query_pos.pl - Will map SNP positions on reference sequence from raw blast results to query postions.

=head1 SYNOPSIS

    my %options;
    my $results = GetOptions (\%options,
			      'query_acc_map|q=s',
			      'raw_blast_list|r=s',
			      'merged_table|m=s',
			      'output_file|o=s'
                              );
 USAGE: map_query_pos.pl
       --query_acc_map=/path/to/some/query_map.txt
       --raw_blast_list=/path/to/blast_raw.list
       --merged_table=/path/to/snp_merged_table.txt
       --output_file=/path/to/output_file.txt
     [ --log=/path/to/file.log
       --debug=3
       --help
     ]

=head1 OPTIONS

B<--query_acc_map,-q>
    Input tab delimited file with 2 columns, first column is the query acc and the second is full path to fasta file

B<--raw_blast_list,-r>
    Path to list of raw blast files. Expects blast results from extract SNP regions as subject and query fastas as query database.

B<--merged_table,-m>
    Path to SNP merged table. Created form SNP-verify or Skirret pipelies.

B<--output_file,-o>
    Path to output file which will be written.

B<--log,-l>
    Logfile.

B<--debug,-d>
    1,2 or 3. Higher values more verbose. [Default:1]

B<--help,-h>
    Print this message

=head1  DESCRIPTION

    This script will read in Blast results and attempt to map positions in some reference to query fasta. Will read in blast
    results and find the best mapping (sometimes a tie). 
 
=head1  INPUT
    
    See SNP::MergedTable for a description of MergedTable format (--merged_table option).

=head1 OUTPUT
    

=head1  CONTACT

    Kevin Galens
    kgalens@gmail.com

=cut
