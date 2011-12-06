package SNP::MergedTable::Row;

use strict;
use warnings;
use Carp;
use Data::Dumper;

our $AUTOLOAD;

my %columns = ( 'molecule'      => { 'index' => 0,  'req' => 1, 'multiple' => 0 }, 
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
				'properties'    => { 'index' => 17, 'req' => 0, 'multiple' => 0 }
);

sub new {
  my ($class, $row) = @_;
  my $self = { 
			  'row' => {},
			  'queries' => [],
			  'query_bases' => {},
			  'num_hits' => {}
			 };
  $self->row( $row ) if( defined( $row ) );
  bless($self, $class);
  return $self;
}

sub row {
  my ($self, $row) = @_;
  $self->{'row'} = $row if( defined( $row ) );
  return $self->{'row'};
}

sub is_column_multiple {
  my ($class, $col) = @_;
  return $columns{$col}->{'multiple'};
}

sub query_base {
  my ($self, $query, $base) = @_;
  croak("Subroutine query_base requires query name") unless( defined( $query ) );
  $self->{'query_bases'}->{$query} = $base if( defined( $base ) );
  return $self->{'query_bases'}->{$query} or undef;
}

sub get_queries {
  my ($self) = @_;
  return keys %{$self->{'query_bases'}};
}

sub num_hit {
  my ($self, $query, $num_hit) = @_;
  croak("Subroutine num_hit requires query name") unless( defined( $query ) );
  $self->{'num_hits'}->{$query} = $num_hit if( defined( $num_hit ) );
  return $self->{'num_hits'}->{$query} or undef;
}

sub to_obj {
  my ($string, @queries) = @_;
 #  my @c = split(/\t/, $string);
#   my @row;
  
#   foreach my $col ( &get_sorted_columns() ) {
# 	if( &is_columns_multiple( $col ) ) {
# 	  # Add columns for multiple
# 	  my @tmp = @c[$columns{$col}->{'index'} .. ($columns{$col}->{'index'} + $num_queries - 1)];
# 	  push(@row, \@tmp);

# 	  ## Slice the rows out so indexes will correctly match
# 	  @c = @c[0..$columns{$col}->{'index'},($columns{$col}->{'index'}+$num_queries)..(scalar(@c))];

# 	} else {
# 	  push(@row, $c[$columns{$col}->{'index'}]);
# 	}
#   }

#   new MergedTable::Row( \@row );
}

sub to_string {
  my ($self, $options) = @_;

  ########## Deal with the options ##########
  foreach my $req ( qw(queries) ) {
	die("Option 'queries' is required to ".__PACKAGE__."::to_string")
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

	if( $self->is_column_multiple( $col ) ) {
	  my $missing_value;
	  if( $col eq 'query_bases' ) {
		$missing_value = ($missing_query_bases_nohit) ? "No Hit" : $self->{'row'}->{'refbase'};
	  } elsif( $col eq 'num_hits' ) {
		$missing_value = 0;
	  } else {
		croak("Multiple column $col is not supported");
	  }
	  push(@retval, $self->multiple_column_to_string( $col, $queries, $missing_value ) );
	} else {
	  push(@retval, $self->{'row'}->{$col});
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

  return ($valid, $message);
  
}

sub AUTOLOAD {
  my $self = shift;
  my $name = $AUTOLOAD;
  $name =~ s/.*://;
  
  return if( $name eq 'DESTROY' );
  croak("Method missing: $AUTOLOAD") unless( exists( $columns{$name} ) );
  
  if( @_ ) {
	if( $columns{$name}->{'multiple'} == 1 ) {
	  my $aref = shift;
	  croak("Subroutine $AUTOLOAD requires one argument, an arrayref.") 
		unless( ref( $aref ) eq 'ARRAY' );
	  $self->{'row'}->{ $name } = $aref;
	} else {
	  return $self->{'row'}->{ $name } = shift;
	}
  } else {
	return $self->{'row'}->{ $name };
  }
  
}

1;
