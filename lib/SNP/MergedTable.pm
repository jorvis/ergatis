package SNP::MergedTable;

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
 properties:      Column of key-value pairs in the format of <key>=<value> (ex. verified=false) separated by
                  semi-colons. Optional.

=head1 Synopsis

 my $mtable = SNP::MergedTable::parse( "/path/to/merged.table" );

 # Need to write this part yet.
 $mtable->filter( <STUFF HERE> );

 $mtable->print_to_fh( $fh );

 or

 my $mtable = new SNP::MergedTable();
 $mtable->add_row( new SNP::MergedTable::Row( $row ) );
 $mtable->print_to_file( "/path/to/output/merged.table" );

=head1 Author

Kevin Galens (kgalens@gmail.com)

=cut

use strict;
use warnings;
use Carp;
use SNP::MergedTable::Row;
use SNP::MergedTable::Filter;
use Data::Dumper;

sub new {
  my ($class, %options) = @_;
  my $self = { 
			  'rows' => {},
			  'queries' => [],
			  'query_hash' => {},
			  'total_rows' => 0,
			  'missing_query_bases_nohit' => 0,
			  'include_num_hits' => 0
			 };

  $self->{'missing_query_bases_nohit'} = $options{'missing_query_bases_nohit'} 
	if( exists( $options{'missing_query_bases_nohit'} ) );
  $self->{'include_num_hits'} = $options{'include_num_hits'} 
	if( exists( $options{'include_num_hits'} ) );

  bless( $self, $class );
  return $self;
}

=head1 Methods

=head2 Class Methods

=over 2

=head3 SNP::MergedTable::parse( $file_path )

=over 3

B<Description:> Creates a MergedTable object from a file path.

B<Parameters:> $file_path: path to the merged table file to parse

B<Returns:> a SNP::MergedTable object reference

=back

=cut 

sub parse {
  my ($file) = @_;

  die("Merged table file [$file] does not exist") unless( -e $file );

  my $self = new SNP::MergedTable();

  # Parse the header to determine the queries.
  open(IN, "< $file") or die("Can't open $file: $!");  
  chomp( my $header_line = <IN> );
  my @header = split(/\t/, $header_line);
  
  # Grab the index for the query_bases
  my $query_bases_index = SNP::MergedTable::Row::get_index_by_column_name( 'query_bases' );
  croak("Could not get column index for header value 'query_bases") unless( defined( $query_bases_index ) );

  # Grab the next header column name (so we know what to look for)
  my $next_column_name = SNP::MergedTable::Row::get_column_name_by_index( $query_bases_index + 1 );
  croak("Could not grab next column name with index ".($query_bases_index + 1) ) unless( defined( $next_column_name ) );

  # Finally grab the query names
  foreach my $head ( @header[$query_bases_index .. scalar(@header)] ) {
	last if( $head eq $next_column_name );
	$self->add_query( $head );
  }

  # Now we can start parsing the rows
  my $num_hits_flag = 0;
  while( my $row_string = <IN> ) {
	chomp $row_string;
	my $row = SNP::MergedTable::Row::to_obj( $row_string, $self->queries() );
	unless( $num_hits_flag ) {
	    $self->include_num_hits(0);
	    my %h = $row->get_all_num_hits();
	    $self->include_num_hits(1) if( keys %h );
	    $num_hits_flag = 1;
	}
	$self->add_row( $row );
  }
  
  return $self;
}

=back

=head2 Instance Methods

=head3 $merged_table->missing_query_bases_nohit( $var )

=over 3

B<Description:> Sets whether a missing query_base should have a value of 'No Hit'. By
 default will use the reference base if missing. If no parameter is given, will return
 current value set for the parameter.

B<Parameters:> $var: A non-zero value indicates 'No Hit' value should be used for missing query bases.

B<Returns:> The value for the missing_query_bases_nohit param

=back

=cut 
sub missing_query_bases_nohit {
  my ($self, $var) = @_;
  $self->{'missing_query_bases_nohit'} = $var if( defined( $var ) );
  return $self->{'missing_query_bases_nohit'};
}

=head3 $merged_table->include_num_hits( $var )

=over 3

B<Description:> If set to non-zero, a num_hit column for each query will be included in 
 the output file. If set to zero, these columns will be omitted.

B<Parameters:> $var: A non-zero value indicates num_hit columns should be included in the output.

B<Returns:> The value for the include_num_hits parameter.

=back

=cut 
sub include_num_hits {
  my ($self, $var) = @_;
  $self->{'include_num_hits'} = $var if( defined( $var ) );
  return $self->{'include_num_hits'};
}

=head3 $merged_table->queries( [$queries] )

=over 3

B<Description:> Used to set the queries. If an array ref of queries is provided,
 all queries previously set will be removed and the new set added. If no parameters 
 are passed in, will return the current list of queries.

B<Parameters:> $queries: used to set the query list.

B<Returns:> Array of query names

=back

=cut 
sub queries {
  my ($self, $queries) = @_;
  
  if( defined( $queries ) ) {
	die("Parameter queries should be an array reference") unless( ref( $queries ) eq 'ARRAY' );

	# If this was called, will reset all the queries, so remove them all first
	$self->_remove_queries();

	# We are calling add_query here instead of just storing the array reference
	# because we also need to count and store the queries in a hash. Keep the storage
	# of queries consistent (even if not optimal).
	foreach my $q ( @{$queries} ) {
	  $self->add_query( $q );
	}
  }
  return @{$self->{'queries'}};
}

sub _remove_queries {
  my ($self) = @_;
  $self->{'queries'} = [];
  $self->{'query_hash'} = {};
}

=head3 $merged_table->add_query( $query )

=over 3

B<Description:> Will add a query to this table.

B<Parameters:> $query: the name of the query.

B<Returns:> 1 on success

=back

=cut 
sub add_query {
  my ($self, $query) = @_;
  
  # Push is on the array to preserve the order
  push(@{$self->{'queries'}}, $query);

  # Also store it in the hash for quicker retrieval
  $self->{'query_hash'}->{$query} = 0 unless( exists( $self->{'query_hash'}->{$query} ) );
  $self->{'query_hash'}->{$query}++;

  return 1;
}

=head3 $merged_table->num_queries( )

=over 3

B<Description:> Returns the total number of queries added to this table.

B<Parameters:> None

B<Returns:> Number of queries.

=back

=cut 
sub num_queries {
  my ($self) = @_;
  return scalar( @{$self->{'queries'}} );
}

=head3 $merged_table->query_exists( $query_name )

=over 3

B<Description:> Returns true if the query exists in this table.

B<Parameters:> $query: query name

B<Returns:> Will return the number of times this query has been added to the file.

=back

=cut 
sub query_exists {
  my ($self, $query) = @_;
  my $retval = exists( $self->{'query_hash'}->{$query} ) ? $self->{'query_hash'}->{$query} : 0;
  return $retval;
}

=head3 $merged_table->add_row( $row, $queries )

=over 3

B<Description:> Will accept a string representation of a row or a SNP::MergedTable::Row object
 and store it as a row in the table. 

 *NOTE: If a string representation of a row is passed in, the current list of queries is
 used to determine the column headers unless specified with the $queries parameter. 
 If a row string is passed in before the queries are set and the $queries parameter is not set,
 this method will fail.

B<Parameters:> $row: Either a string representation of a row or a SNP::MergedTable::Row object reference
 $queries: array reference of query names used when $row is a string. Indicates the headers for the query_bases
  and num_hits columns.

B<Returns:> Nothing

=back

=cut 
sub add_row {
    my ($self, $row, $queries) = @_;
    
    my $obj = $row;
    if( ref( $row ) ne 'SNP::MergedTable::Row' ) {
        $queries = $self->queries() unless( defined( $queries ) );
        $obj = SNP::MergedTable::Row::to_obj( $queries );
    }

    $self->{'rows'}->{$row->molecule} = {} unless( exists( $self->{'rows'}->{$row->molecule} ) );
    $self->{'rows'}->{$row->molecule}->{$row->refpos} = {}
        unless( exists( $self->{'rows'}->{$row->molecule}->{$row->refpos} ) );
    $self->{'rows'}->{$row->molecule}->{$row->refpos}->{$row->gene_name} = $row;
    
    $self->{'total_rows'}++;

}

=head3 $merged_table->remove_row( )

=over 3

B<Description:> Will remove a row (or rows) if they exist. If row specified does
    not exist, program will die. 

B<Parameters:> $row: a SNP::MergedTable::Row object.

B<Returns:> The number of rows removed from the table (should only ever be 1).

=back

=cut
sub remove_row {
    my ($self, $row) = @_;
    die("remove_row requires a SNP::MergedTable::Row object")
        unless( defined( $row ) && ref( $row ) eq 'SNP::MergedTable::Row' );
        
    my $gene;
    $gene = $row->gene_name unless( $row->gene_name eq 'intergenic' );
    my $count = $self->remove_row_by_position( $row->molecule, $row->refpos, $gene );
    die("Expected to remove 1 and only 1 row. Removed $count rows.") unless( $count == 1 );
    
    return $count;
}

=head3 $merged_table->remove_row_by_position( )

=over 3

B<Description:> Will remove a row (or rows) if they exist. If row specified does
    not exist, program will die. If a gene name is provided, only the SNP located
    on that specified gene will be removed. If no gene name is provided, multiple
    rows might be removed if multiple SNPs occur at the specified refpos.

B<Parameters:> $molecule: Molecule name on which the SNP is located.
               $refpos: The position on the molecule.
               $gene: [Optional] The gene name on which the SNP is located.

B<Returns:> The number of rows removed from the table.

=back

=cut
sub remove_row_by_position {
    my ($self, $molecule, $refpos, $gene) = @_;
    die("Cannot remove row [$molecule, $refpos, $gene]. Does not exist")
        unless( $self->row_exists( $molecule, $refpos, $gene ) );
        
    my $count = 0;
    if( defined( $gene ) ) {
        delete( $self->{'rows'}->{$molecule}->{$refpos}->{$gene} );
        $count++;
    } else {
        my @genes = keys %{$self->{'rows'}->{$molecule}->{$refpos}};
        delete( $self->{'rows'}->{$molecule}->{$refpos} );
        $count = scalar(@genes);
    }
    
    $self->{'total_rows'} -= $count;
    return $count;
}

=head3 $merged_table->row_exists( )

=over 3

B<Description:> Checks to see if a specified row exists. If gene name is provided
    will specifically look for a SNP stored under that gene name. If no gene name 
    is provided, will return non-zero if any SNP is present at that refpos.

B<Parameters:> $molecule: Molecule name on which the SNP is located.
               $refpos: The position on the molecule.
               $gene: [Optional] The gene name on which the SNP is located.

B<Returns:> Non-zero if a row exists describing the SNP.

=back

=cut
sub row_exists {
    my ($self, $molecule, $refpos, $gene) = @_;
    my $retval = 0;
    
    if( exists( $self->{'rows'}->{$molecule}->{$refpos} ) ) {
        if( defined( $gene ) ) {
            $retval = 1 if( exists( $self->{'rows'}->{$molecule}->{$refpos}->{$gene} ) );
        } else {
            $retval = 1;
        }
    }
    return $retval;
}

=head3 $merged_table->get_rows( )

=over 3

B<Description:> Retrieves the rows for the table.

B<Parameters:> None

B<Returns:> An array of rows, sorted by molecule name then reference SNP position

=back

=cut 
sub get_rows {
  my ($self) = @_;
  my @retval;
  
  foreach my $mol ( sort { $a cmp $b } keys %{$self->{'rows'}} ) {
      foreach my $refpos ( sort { $a <=> $b } keys %{$self->{'rows'}->{$mol}} ) {
          foreach my $gene_name ( sort { $a cmp $b } keys %{$self->{'rows'}->{$mol}->{$refpos}} ) {
              push(@retval, $self->{'rows'}->{$mol}->{$refpos}->{$gene_name});
          }
      }
  }
  
  return @retval;
}

=head3 $merged_table->get_row_by_position( $molecule, $refpos )

=over 3

B<Description:> Retrieve a row by the SNP position on a reference molecule

B<Parameters:> $molecule: name of the molecule on which the snp is located
 $refpos: position on that molecule.

B<Returns:> Returns a SNP::MergedTable::Row object reference is row exists for a SNP
 at that position, otherwise will return undef.

=back

=cut 
sub get_row_by_position {
  my ($self, $molecule, $refpos) = @_;
  my @retval;
  if( exists( $self->{'rows'}->{$molecule} ) && exists(  $self->{'rows'}->{$molecule}->{$refpos} ) ) {
      @retval = values( %{$self->{'rows'}->{$molecule}->{$refpos}} );
  }
  return @retval;
}

=head3 $merged_table->get_header( )

=over 3

B<Description:> Retrieve header text. Will die if no queries have been set.

B<Parameters:> None

B<Returns:> An array containing the names for the columns

=back

=cut 
sub get_header {
  my ($self) = @_;

  my @header;
  foreach my $col ( SNP::MergedTable::Row->get_sorted_columns() ) {
	next if( $col eq 'num_hits' && !$self->include_num_hits );
	if( SNP::MergedTable::Row::is_column_multiple( $col ) ) {
	  $self->_determine_queries() unless( $self->num_queries > 0 );
	  my @queries = $self->queries;
	  if( $col eq 'num_hits' ) {
	      map { $_ = "num_hits: $_" } @queries;
	  }
	  push( @header, @queries );
	} else {
	  push( @header, $col );
	}
	
  }
  @header;
}

sub _determine_queries {
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

=head3 $merged_table->filter( $fh )

=over 3

B<Description:> Will filter the merged table file. Will only work with current
    rows. Any future rows added will not be filtered.

B<Parameters:> $filter_name: the name of the filter. Must be valid

B<Returns:> Nothing.

=back

=cut
sub filter {
    my ($self, $filter_name, @args) = @_;    
    SNP::MergedTable::Filter::filter( $filter_name, $self, @args );
}

=head3 $merged_table->print_to_fh( $fh )

=over 3

B<Description:> Will write Merged Table to the passed in file handle

B<Parameters:> $fh: File handle

B<Returns:> None

=back

=cut 
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

=head3 $merged_table->print_to_file( $file_path )

=over 3

B<Description:> Will write Merged Table to the passed in file

B<Parameters:> $file_path: path to output file

B<Returns:> None

=back

=cut 
sub print_to_file {
  my ($self, $file) = @_;
  open( my $fh, "> $file") or croak("Could not open $file for writing: $!");
  $self->print_to_fh( $fh );
  close( $fh );
  return 1;
}

1;
