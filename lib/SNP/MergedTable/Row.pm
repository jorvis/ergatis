package SNP::MergedTable::Row;

use strict;
use warnings;
use Carp;
use Data::Dumper;

our $AUTOLOAD;

my %columns = ( 
		'molecule'      => { 'index' => 0,  'req' => 1, 'multiple' => 0 }, 
		'refpos'        => { 'index' => 1,  'req' => 1, 'multiple' => 0 },
		'syn'           => { 'index' => 2,  'req' => 1, 'multiple' => 0 },
		'refbase'       => { 'index' => 3,  'req' => 1, 'multiple' => 0 },
		'query_bases'   => { 'index' => 4,  'req' => 0, 'multiple' => 1 },
		'gene_name'     => { 'index' => 5,  'req' => 1, 'multiple' => 0 }, 
		'product'       => { 'index' => 6,  'req' => 1, 'multiple' => 0 },
		'gene_start'    => { 'index' => 7,  'req' => 1, 'multiple' => 0 },
		'gene_stop'     => { 'index' => 8,  'req' => 1, 'multiple' => 0 },
		'gene_length'   => { 'index' => 9,  'req' => 1, 'multiple' => 0 },
		'snps_per_gene' => { 'index' => 10, 'req' => 1, 'multiple' => 0 },
		'pos_in_gene'   => { 'index' => 11, 'req' => 1, 'multiple' => 0 },
		'ref_codon'     => { 'index' => 12, 'req' => 1, 'multiple' => 0 },
		'ref_aa'        => { 'index' => 13, 'req' => 1, 'multiple' => 0 },
		'snp_codon'     => { 'index' => 14, 'req' => 1, 'multiple' => 0 },
		'snp_aa'        => { 'index' => 15, 'req' => 1, 'multiple' => 0 },
		'num_hits'      => { 'index' => 16, 'req' => 0, 'multiple' => 1 },
		'properties'    => { 'index' => 17, 'req' => 0, 'multiple' => 0 },
		'other'         => { 'index' => 18, 'req' => 0, 'multiple' => 0 },
		);

sub new {
  my ($class, $row) = @_;
  my $self = { 
			  'row' => {},
			  'queries' => [],
			  'query_bases' => {},
			  'num_hits' => {}
			 };
  bless($self, $class);
  $self->row( $row ) if( defined( $row ) );
  return $self;
}

sub row {
  my ($self, $row) = @_;
  $self->{'row'} = $row if( defined( $row ) );
  return $self->{'row'};
}

sub is_column_multiple {
  my ($col) = @_;
  die("Could not find col: $col") unless( exists( $columns{$col} ) );
  return $columns{$col}->{'multiple'};
}

sub query_base {
  my ($self, $query, $base) = @_;
  croak("Subroutine query_base requires query name") unless( defined( $query ) );
  $self->{'query_bases'}->{$query} = $base if( defined( $base ) );
  return $self->{'query_bases'}->{$query} or undef;
}

sub num_hit {
  my ($self, $query, $num_hit) = @_;
  croak("Subroutine num_hit requires query name") unless( defined( $query ) );
  $self->{'num_hits'}->{$query} = $num_hit if( defined( $num_hit ) );
  return $self->{'num_hits'}->{$query};
}

sub get_all_num_hits {
    my ($self) = @_;
    return %{$self->{'num_hits'}};
}

sub get_queries {
  my ($self) = @_;
  return keys %{$self->{'query_bases'}};
}

sub to_obj {
  my ($string, @queries) = @_;

  # Delimited by tabs
  my @c = split(/\t/, $string);

  # Count the number of queries
  my $num_queries = scalar(@queries);

  # We should make sure that we have at least the correct number of columns
  # based on the queries passed in. 
  my @all_columns = &get_columns();
  my $num = (@all_columns-2) + $num_queries;
  
  # properties and other columns can be missing.
  $num -= 2;

  die("Row [$c[0], $c[1], $c[2]] did not have enough columns. Expected at least $num but found: ".scalar(@c).
      ". Number of queries: $num_queries")
      if( $num > @c );
      
  # Create the object we will eventually return
  my $row = new SNP::MergedTable::Row;

  # Grab all the columns and add the values to the object
  foreach my $col ( &get_sorted_columns() ) {

      if( &is_column_multiple( $col ) ) {

 	  # Add columns for multiple
 	  my @tmp = @c[$columns{$col}->{'index'} .. ($columns{$col}->{'index'} + $num_queries - 1)];
	  die("Didn't slice correctly") unless( @tmp == $num_queries );
	  for( my $i = 0; $i < scalar(@tmp); $i++ ) {
	      if( $col eq 'query_bases' ) {
		  $row->query_base( $queries[$i], $tmp[$i] );
	      } elsif( $col eq 'num_hits' ) {
		  $row->num_hit($queries[$i], $tmp[$i]);
	      }
	  }

 	  ## Slice the rows out so indexes will correctly match
 	  @c = @c[0..$columns{$col}->{'index'},($columns{$col}->{'index'}+$num_queries)..(scalar(@c))];
      } elsif( $col eq 'other' ) {
	  # This is the last column so slurp everything up
	  my $len = scalar(@c);
	  my @rest = grep { defined } @c[$columns{$col}->{'index'}..($len-1)];
	  $row->set_column( $col, join("\t", @rest) );
      } else {
	  my $i = $columns{$col}->{'index'};
	  if( $i < @c && defined( $c[$i] ) ) {
	      $row->set_column( $col, $c[$columns{$col}->{'index'}] );
	  } elsif( $columns{$col}->{'req'} ) {
	      die("Required column $col was missing");
	  }
      }
  }

  return $row;

}

sub to_string {
  my ($self, $options) = @_;

  ########## Deal with the options ##########
  foreach my $req ( qw(queries) ) {
	die("Option '$req' is required to ".__PACKAGE__."::to_string")
	  unless( exists( $options->{$req} ) );
  }

  my $queries = $options->{'queries'};
  
  my $missing_query_bases_nohit = 0;
  $missing_query_bases_nohit = $options->{'missing_query_bases_nohit'} if( exists( $options->{'missing_query_bases_nohit'} ) );

  my $include_num_hits = 0;
  $include_num_hits = $options->{'include_num_hits'} if( exists( $options->{'include_num_hits'} ) );
  ###########################################
  
  # Make sure we have a valid row before outputting the string
  my ($valid, $msg) = $self->is_row_valid();
  croak("Row is not valid: $msg") unless( $valid );

  my @retval;
  foreach my $col ( &get_sorted_columns() ) {
      next if( $col eq 'num_hits' && !$include_num_hits );

      if( &is_column_multiple( $col ) ) {
	  my $missing_value;
	  if( $col eq 'query_bases' ) {
	      $missing_value = ($missing_query_bases_nohit) ? "No Hit" : $self->{'row'}->{'refbase'};
	  } elsif( $col eq 'num_hits' ) {
	      $missing_value = 0;
	  } else {
	      croak("Multiple column $col is not supported");
	  }
	  my $val = $self->multiple_column_to_string( $col, $queries, $missing_value ) || "";
	  push(@retval, $val);
      } else {
	  my $val = $self->{'row'}->{$col} || "";
	  push(@retval, $val);
      }
  }

  join("\t", @retval);
}

sub multiple_column_to_string {
  my ($self, $col, $queries, $missing_value) = @_;
  my @retval;

  foreach my $query ( @{$queries} ) {
	my $tmp = $self->{$col}->{$query};
	$tmp = $missing_value if ( !defined( $tmp ) );
	push(@retval, $tmp);
  }

  return join("\t", @retval);
}

sub get_sorted_columns {
  my ($self) = @_;
  sort { $columns{$a}->{'index'} <=> $columns{$b}->{'index'} } &get_columns();
}

sub get_columns {
  my ($self) = @_;
  return keys %columns;
}

sub get_index_by_column_name {
  my ($column_name) = @_;
  return unless( exists( $columns{$column_name} ) );
  return $columns{$column_name}->{'index'};
}

sub get_column_name_by_index {
  my ($index) = @_;

  my $retval;
  foreach my $column ( &get_columns() ) {
	if( $columns{$column}->{'index'} == $index ) {
	  $retval = $column;
	  last;
	}
  }
  return $retval;
  
}

## Checks to see that the row is valid. 
sub is_row_valid {
  my ($self) = @_;
  my $valid = 1;
  my $message;

  my $row = $self->{'row'};
  foreach my $col ( keys %columns ) {

	# If it's required but not defined or wasn't provided.
	if( $columns{ $col }->{'req'} && (!exists( $row->{ $col } ) || !defined( $row->{ $col } ) ) ) {
	  $valid = 0;
	  $message = "Column $col is required by not provided";
	} elsif( !$columns{ $col }->{'req'} && !defined( $row->{ $col } ) ) {
	  $row->{ $col } = '';
	}

	last unless( $valid );
  }
  
  # Check to make sure synonomous and non-synonomous labels make sense.
  if( $valid ) {
      unless( $self->gene_name eq 'intergenic' ) {

	  ## We can have multiple SYN,NSYN definitions because of the following two reasons:
	  ## 1. A query aligned multiple times at this reference position
	  ## 2. There are 2 genes in different frames on the reference. This would cause 2 different
	  ##    codons and the same SNP could cause a SYN mutation in on of the genes and a NSYN in 
	  ##    the other.
	  ## If we have only 1 ref_aa, it's the first case and if we have multiple ref_aa, then it's the second.
	  ## I'm not really checking for a combination of the two, which I probably should be.
          my $syn = $self->syn();
          my @syn = split(/\//, $syn );
          my $ref_aa = $self->ref_aa();
	  my @ref_aa = split(/\//, $ref_aa);
          my $snp_aa = $self->snp_aa();
          my @snpaa = split(/\//, $snp_aa );

	  my $raa;
          for( my $i = 0; $i < @syn; $i++ ) {
              die("SNPAA undefined") if( !defined( $snpaa[$i] ) );
	      
	      $raa = $ref_aa[$i] if( $i < scalar(@ref_aa) );
	      
	      if( $raa ne $snpaa[$i] && $syn[$i] eq 'SYN' ) {
		  $syn[$i] = 'NSYN';
	      } elsif( $raa eq $snpaa[$i] && $syn[$i] eq 'NSYN' ) {
		  $syn[$i] = 'SYN';
	      }
              #if( ( $syn[$i] eq 'SYN' && $raa ne $snpaa[$i] ) || ( $syn[$i] eq 'NSYN' && $raa eq $snpaa[$i] )) {
              #    $valid = 0;
              #    my $pos = $self->refpos();
              #    $message = "SYN/NSYN not correct. [Position: $pos, $syn[$i], RefAA: $raa, SNPAA: $snpaa[$i]]";
              #    #print Dumper( $self );
              #}
          }
	  my $new_syn_string = join("/", @syn);
	  if( $new_syn_string ne $syn ) {
	      $self->syn($new_syn_string);
	  }
      }
  }

  return ($valid, $message);
  
}

sub set_column {
  my ($self, $name, $value) = @_;
  
  if( $columns{$name}->{'multiple'} == 1 ) {
	my $aref = $value;
	croak("Subroutine $AUTOLOAD requires one argument, an arrayref.") 
	  unless( ref( $aref ) eq 'ARRAY' );
	$self->{'row'}->{ $name } = $aref;
  } else {
	$self->{'row'}->{ $name } = $value;
  }
  return $self->{'row'}->{ $name };
}

sub AUTOLOAD {
  my $self = shift;
  my $name = $AUTOLOAD;
  $name =~ s/.*://;
  
  return if( $name eq 'DESTROY' );
  croak("Method missing: $AUTOLOAD") unless( exists( $columns{$name} ) );
  
  if( @_ ) {
	my $value = shift;
	$self->set_column( $name, $value );
  } else {
	return $self->{'row'}->{ $name };
  }
  
}

1;
