package SNP::MergedTable;

use strict;
use warnings;
use Carp;
use SNP::MergedTable::Row;

sub new {
  my ($class, %options) = @_;
  my $self = { 
			  'snps' => {},
			  'queries' => {},
			  'missing_query_bases_nohit' => 0,
			  'include_num_hits' => 0
			 };

  bless( $self, $class );
  return $self;
}

sub parse {
  my ($self, $file) = @_;
  warn "sub parse not written yet\n";
  return;

  open(IN, "< $file") or die("Can't open $file: $!");
  while ( my $line = <IN> ) {
	chomp($line);
  }
  close(IN);
}

sub missing_query_bases_nohit {
  my ($self, $var) = @_;
  $self->{'missing_query_bases_nohit'} = $var if( defined( $var ) );
  return $self->{'missing_query_bases_nohit'};
}

sub include_num_hits {
  my ($self, $var) = @_;
  $self->{'include_num_hits'} = $var if( defined( $var ) );
  return $self->{'include_num_hits'};
}

sub queries {
  my ($self) = @_;
  return keys %{$self->{'queries'}};
}

sub add_query {
  my ($self, $query) = @_;
  $self->{'queries'}->{$query} = 1;
}

sub num_queries {
  my ($self) = @_;
  return scalar( keys %{$self->{'queries'}} );
}

sub query_exists {
  my ($self, $query) = @_;
  return exists( $self->{'queries'}->{$query} );
}

sub determine_queries {
  my ($self) = @_;
  my @rows = $self->get_rows;
  foreach my $r ( @rows ) {
	foreach my $q ( $r->get_queries ) {
	  unless( $self->query_exists( $q ) ) {
		$self->add_query( $q );
	  }
	}
  }
}

sub add_row {
  my ($self, $row) = @_;
  $self->{'snps'}->{$row->molecule} = {} unless( exists( $self->{'snps'}->{$row->molecule} ) );
  $self->{'snps'}->{$row->molecule}->{$row->refpos} = $row;
}

sub get_row_by_position {
  my ($self, $molecule, $refpos) = @_;
  my $retval;
  if( exists( $self->{'snps'}->{$molecule} ) && exists(  $self->{'snps'}->{$molecule}->{$refpos} ) ) {
	$retval = $self->{'snps'}->{$molecule}->{$refpos};
  }
  return $retval;
}

sub get_header {
  my ($self) = @_;

  ## Before we get the header, make sure we've collected all the queries
  ## present across all rows.
  $self->determine_queries();

  my @header;
  foreach my $col ( SNP::MergedTable::Row->get_sorted_columns() ) {
	next if( $col eq 'num_hits' && !$self->include_num_hits );
	if( SNP::MergedTable::Row->is_column_multiple( $col ) ) {
	  die("Could not find any queries in any of the added rows") unless( $self->num_queries > 0 );
	  push( @header, $self->queries );
	} else {
	  push( @header, $col );
	}
	
  }
  @header;
}

sub get_rows {
  my ($self) = @_;
  my @retval;
  foreach my $mol ( sort { $a cmp $b } keys %{$self->{'snps'}} ) {
	foreach my $refpos ( sort { $a <=> $b } keys %{$self->{'snps'}->{$mol}} ) {
	  push(@retval, $self->{'snps'}->{$mol}->{$refpos});
	}
  }
  return @retval;
}

sub print_to_fh {
  my ($self, $fh) = @_;

  ## Print the header.
  my @header = $self->get_header();
  print $fh join("\t", @header);
  print $fh "\n";
  
  my @queries = $self->queries;
  my $options = { 'queries' => \@queries,
				  'missing_query_bases_nohit' => $self->missing_query_bases_nohit(),
				  'include_num_hits' => $self->include_num_hits()
				};
				
  foreach my $row ( $self->get_rows ) {
	print $fh $row->to_string( $options )."\n";
  }
  return 1;
}

sub print_to_file {
  my ($self, $file) = @_;
  open( my $fh, "> $file") or croak("Could not open $file for writing: $!");
  $self->print_to_fh( $fh );
  close( $fh );
  return 1;
}

1;

=head1 Module

SNP::MergedTable.pm - Module for writing/parsing SNP MergedTable files

=head1 Description


This module is used for writing/parsing a format used with SNP descriptions, MergedTable. The files is tab delimited and has the following columns:

 molecule:        An identifier for the molecule.
 refpos:          Position of SNP in reference
 syn?:            SYN or NSYN if within gene, 'NA' otherwise
 refbase:         The base in the reference
 <queries>:       Variable number of columns, one for each query. The value will be the base of the SNP in the query. 
                  If the region was not found in the query, the string 'No Hit' will be here.
 gene_name:       If the SNP is with a gene, the name of that gene, otherwise the string 'None'
 product:         The common name of the gene this SNP is within. The string "intergenic" if not in a gene
 gene_start:      Start of the gene.
 gene_stop:       Stop of the gene.
 gene_length:     Length of gene if SNP is within gene.
 snps_per_gene:   The number of SNPs in this gene.
 pos_in_gene:     Position of SNP within gene.
 ref_codon:       The reference codon
 ref_aa:          The reference amino acid.
 query_codon:     The query codon (if multiple, separated by /, example ACT/AGT)
 query_aa:        The query amino acid (if multiple, separatedy by /, example T/S)
 <num_hits>:      A column for each query, listing the number of hits from blast results [used in snp-verify]. Optional.
 properties:      Column of key-value pairs in the format of <key>=<value> (ex. verfified=false) separated by
                  semi-colons. Optional.

=head1 Author

Kevin Galens (kgalens@gmail.com)

=cut
