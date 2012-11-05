#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Pod::Usage;
use Data::Dumper;
use SNP::MergedTable;
use SNP::MergedTable::Row;

our $MISSING_QUERIES_NO_HIT = 0;

################# GLOBALS AND CONSTANTS ###################
my $debug = 1;
my ($ERROR, $WARN, $DEBUG) = (1,2,3);
my $logfh;
my $snps_per_gene = {};
my $snp_nucleotides = {};
my $map = 0;
my $query_map = {};
############################################################

my %options; 
my $results = GetOptions(\%options,
			 "input_file|i=s",
			 "output_file|o=s",
			 "query_map|q=s",
			 "log|l=s",
			 "debug|d=s",
			 "help|h"
			 );

&check_options(\%options);

print "Parsing nucmer snp file: $options{'input_file'}\n";
my $table = &parse_nucmer_snp_file( $options{'input_file'} );
print "Done parsing input file. Printing to output file: $options{'output_file'}\n";
$table->print_to_file( $options{'output_file'} );
print "Finished. Output here: $options{'output_file'}\n";

sub parse_nucmer_snp_file {
  my ($input_file) = @_;

  my $table = new SNP::MergedTable;

  open(my $ifh, "< $input_file") or &_log($ERROR, "Could not open $input_file: $!");
  while ( my $line = <$ifh> ) {
	chomp $line;
	next if( $line =~ /p1/ );
	&_log($ERROR, "Expected line to start with SNP position [$line]") 
	  unless( $line =~ /^\d+/ );
	my @cols = split(/\t/, $line);
	&create_merged_table_row( $table, @cols );
  }
  close($ifh);

  &add_snp_counts( $table );
  return $table;
}

## We add the SNP codon, SNP aa and SYN/NSYN at the very end, in case we have
# multiple variants. This is also the time we count the number of 
# snps per gene.
sub add_snp_counts {
  my ($table) = @_;
  foreach my $row ( $table->get_rows ) {

	## Add the query codon/aa
	unless( $row->gene_name eq 'intergenic' ) {
	  die("Could not find snp codon/aa info for ".$row->refpos) 
		unless( exists( $snp_nucleotides->{$row->refpos} ) );

	  my (@codons, @aas, @syn);
	  foreach my $refbase ( keys %{$snp_nucleotides->{$row->refpos}} ) {
	      foreach my $gene_name ( keys %{$snp_nucleotides->{$row->refpos}->{$refbase}} ) {
	          next unless( $gene_name eq $row->gene_name );
	          my $codon = $snp_nucleotides->{$row->refpos}->{$refbase}->{$gene_name}->{'codon'};
	          my $snpaa = $snp_nucleotides->{$row->refpos}->{$refbase}->{$gene_name}->{'aa'};
	          push(@codons, $codon);
	          push(@aas, $snpaa);
	          my $t = ( $row->ref_aa() eq $snpaa ) ? 'SYN' : 'NSYN';
	          push(@syn, $t);
	      }
	  }
	  $row->snp_codon( join("/", @codons) );
	  $row->snp_aa( join("/", @aas) );
	  $row->syn( join("/", @syn) );
	}

	my @genes = split("/", $row->gene_name);
	my @snps_per_gene = ();
	foreach my $g ( @genes ) {

	  ## Skip the snps in intergenic regions. We don't need to count
	  ## how many snps lie within intergenic regions.
	  if( $g eq 'intergenic' ) {
		push(@snps_per_gene, "NA");
	  } else {
		die("Could not find snps_per_gene for ".$row->molecule." and $g")
		  unless( exists( $snps_per_gene->{$row->molecule}->{$g} ) );
		push(@snps_per_gene, $snps_per_gene->{$row->molecule}->{$g});
	  }

	}
	$row->snps_per_gene( join("/", @snps_per_gene) );
  }
}

sub create_merged_table_row {
  my ($table, @cols) = @_;

  # We might have multiple rows at a certain position if the SNP
  # is found within overlapping genes. This would result in the reference
  # position having multiple codons/amino acids, each being represented by
  # an individual row.
  my @rows = $table->get_row_by_position( $cols[10], $cols[0] );

  # If we passed in a map, replace the query id
  my $query = $cols[11];
  if( $map && exists( $query_map->{$query} ) ) {
      $query = $query_map->{$query};
  }

  ## This is the first time we've seen this position (no rows were returned).
  if ( @rows == 0 ) {

	# Check to see if there is an indel here.
	# merged table doesn't support indels.
	if ( $cols[1] eq '.' || $cols[2] eq '.' ) {
	  die("Input file should not have indels. Please use input file without indels");
	}
	
	# Does this position overlap multiple genes? We should have
	# multiple gene names, gene starts, gene stops, position_in_gene, 
	# product, gene_direction, reference codons and reference amino acids
	# if this is the case.
	my @genes = &_get_multi_values( $cols[12] );
	my @strands = &_get_multi_values( $cols[18] );
	my @lefts = &_get_multi_values( $cols[13] );
	my @rights = &_get_multi_values( $cols[14] );
	my @positions_in_genes = &_get_multi_values( $cols[15] );
	my @products = &_get_multi_values( $cols[17] );
	my @ref_codons = &_get_multi_values( $cols[19] );
	my @ref_aas = &_get_multi_values( $cols[20] );
	my @query_codons = &_get_multi_values( $cols[21] );
	my @query_aas = &_get_multi_values( $cols[22] );
	
	for( my $i = 0; $i < scalar( @genes ); $i ++ ) {
	    my $row = new SNP::MergedTable::Row;

	    $row->molecule( $cols[10] );
	    $row->refpos( $cols[0] );
	    $row->refbase( $cols[1] );
	    $row->query_base( $query, $cols[2] );
	    $row->gene_name( $genes[$i] );
	    
	    if( $row->gene_name eq 'intergenic' ) {
		$row->gene_length( "NA" );
		$row->snp_codon( "NA" );
		$row->snp_aa( "NA" );
		$row->syn("NA");
		$row->gene_start( "NA" );
		$row->gene_stop( "NA" );
		$row->pos_in_gene( "NA" );
		$row->ref_codon( "NA" );
		$row->ref_aa( "NA" );
		$row->product( "NA" );
	    } else {
		
		# For some reason, if we have multiple genes, sometimes multiple strands
		# aren't specified. So just use the first one (maybe this is bad?)
		# Needs to work.
		my $length = $rights[$i] - $lefts[$i] + 1;
		my $strand = ( defined($strands[$i]) ) ? $strands[$i] : $strands[0];
		my ($start, $stop) = ($lefts[$i], $rights[$i]);
		($stop, $start) = ($start, $stop) if( $strand == -1 );
		
		$row->gene_start( $start );
		$row->gene_stop( $stop );
		$row->gene_length( $length );
        	$row->pos_in_gene( $positions_in_genes[$i] );
        	$row->ref_codon( $ref_codons[$i] );
        	$row->ref_aa( $ref_aas[$i] );
        	$row->product( $products[$i]  );
		
		## If the query base is different from reference
		if( $row->refbase ne $cols[2] ) {
		    $snp_nucleotides->{$row->refpos} = {} unless( exists( $snp_nucleotides->{$row->refpos} ) );
		    $snp_nucleotides->{$row->refpos}->{$cols[2]} = {} unless( exists( $snp_nucleotides->{$row->refpos}->{$cols[2]} ) );
		    $snp_nucleotides->{$row->refpos}->{$cols[2]}->{$row->gene_name} = {'codon' => $query_codons[$i],'aa' => $query_aas[$i]};
		}
		
		# Count snps per gene
		$snps_per_gene->{$row->molecule} = {} unless( exists( $snps_per_gene->{$row->molecule} ) );
		$snps_per_gene->{$row->molecule}->{$row->gene_name} = 0
		    unless( exists( $snps_per_gene->{$row->molecule}->{$row->gene_name} ) );
		$snps_per_gene->{$row->molecule}->{$row->gene_name}++;

	    }
	    
	    $table->add_row( $row );
	}
    } else {
	
	for( my $i = 0; $i < scalar(@rows); $i++ ) {
	    my $row = $rows[$i];
	    my $base = $row->query_base( $query );
	    die("Already have query base for query $query [refpos: $cols[0]]") if( defined( $base ) );
	    $row->query_base( $query, $cols[2] );

	    unless( $row->gene_name eq 'intergenic' ) {
		my @gene_names= &_get_multi_values( $cols[12] );
		my @query_codons = &_get_multi_values( $cols[21] );
		my @query_aas = &_get_multi_values( $cols[22] );
		
		for( my $j = 0; $j < scalar(@query_codons); $j++ ) {
		    next unless( $gene_names[$j] eq $row->gene_name );
		    
		    ## If the query base is different from reference
		    if( $row->refbase ne $cols[2] ) {
          		$snp_nucleotides->{$row->refpos} = {} unless( exists( $snp_nucleotides->{$row->refpos} ) );
          		$snp_nucleotides->{$row->refpos}->{$cols[2]} = {} unless( exists( $snp_nucleotides->{$row->refpos}->{$cols[2]} ) );
          		$snp_nucleotides->{$row->refpos}->{$cols[2]}->{$row->gene_name} = {'codon' => $query_codons[$j],'aa' => $query_aas[$j]};
		    }
		}
	    }
	}
    }
}

sub _get_multi_values {

    my ($value, $expect) = @_;
    my @retval = split(/\//, $value);
    die("Expected to find $expect values, found ".scalar(@retval) ) 
        if( defined( $expect ) && scalar(@retval) != $expect );
    return @retval;
}

sub parse_query_map {
    my ($file) = @_;
    open(IN, "< $file" ) or die("Can't open $file: $!");

    my %retval;
    while( my $line = <IN> ) {
	chomp( $line );
	next if( $line =~ /^\#/ || $line =~ /^\s*$/ );
	my @c = split(/\t/, $line);
	my $org_id = shift( @c );
	foreach my $mol ( @c ) {
	    $retval{$mol} = $org_id;
	}
    }
    close(IN);

    return \%retval;
}

sub check_options {

  my $opts = shift;

  if ( $opts->{'help'} ) {
	&_pod;
  }

  if ( $opts->{'log'} ) {
	open( $logfh, "> $opts->{'log'}") or die("Can't open log file ($!)");
  }

  foreach my $req ( qw(input_file output_file) ) {
	&_log($ERROR, "Option $req is required") unless( $opts->{$req} );
  }

  $debug = $opts->{'debug'} if( $opts->{'debug'} );

  $map = 1 if( $opts->{'query_map'} );
  $query_map = &parse_query_map( $opts->{'query_map'} ) if( $opts->{'query_map'} );
}

sub _log {
  my ($level, $msg) = @_;
  if ( $level <= $debug ) {
	print STDOUT "$msg\n";
  }
  print $logfh "$msg\n" if( defined( $logfh ) );
  exit(1) if( $level == $ERROR );
}

sub _pod {
  pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDERR} );
}


=head1 NAME

nucmer_snp2_merged_table.pl - Will take in a concatenated snp file (output of SNP add gene info)

=head1 SYNOPSIS

 nucmer-snp2merged.pl
       --input_file=/path/to/some/input.snps
       --output=/path/to/merged_input.snps
     [ --log=/path/to/file.log
       --debug=3
       --help
     ]

=head1 OPTIONS

B<--input_file,-i> 
    Output from snp-add-gene-info.pl. See INPUT for description of format

B<--output_file,-o>
    Output file is a SNP Merged Table. See OUTPUT section for description of format.

B<--query_map,-q>
    A file which identifies if there are mulitple molecules per query.
    Each lines represents a query and should have the following tab delimited columns:

    query_id     molecule_id-1     molecule_id-2     ... molecule_id-n

B<--log,-l>
    Logfile.

B<--debug,-d>
    1,2 or 3. Higher values more verbose.

B<--help,-h>
    Print this message

=head1  DESCRIPTION

 This merged table format is an internal output format to view snps across a number of query genomes. This script 
 converts the output from snp-add-gene-info into the merged table format. The snp-add-gene-info output 
 format is the same as the output received from nucmer-show-snps program, with a few additional columns 
 added to the end. See the INPUT section for details on the formats.
 
=head1  INPUT

=head2 snp-add-gene-info output

The first 12 columns are the same as those output from show-snps [see show-snps documentation for 
description of columns]:
p1, ref_base, query_base, p2, buff, dist, len_r, len_q, frm1, frm2, ref_contig, query_contig

With these additional columns:
gene_id, position_in_gene, syn_nonsyn, product, gene_direction, ref_codon, ref_amino_acid, query_codon, query_amino_acid, num_homopolymer

=head1 OUTPUT

=head2 Columns in the output:

 See the perldoc for SNP::MergedTable.pm.

=head1  CONTACT

 Kevin Galens
 kgalens@gmail.com

=cut
