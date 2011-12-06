#!/usr/bin/env perl

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

=head1 DEVEL

/export/svn/ergatis/src/perl/nucmer-snp2merged.pl --input_file /usr/local/projects/PVCHO/kgalens/output_repository/snp-add-gene-info/9031_default/snp-add-gene-info.default.no_indels.snps --output_file /tmp/output.file

=cut

use strict;
use warnings;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Pod::Usage;
use Data::Dumper;
use SNP::MergedTable;
use SNP::MergedTable::Row;

our $MISSING_QUERIES_NO_HIT = 0;

############# GLOBALS AND CONSTANTS ################
my $debug = 1;
my ($ERROR, $WARN, $DEBUG) = (1,2,3);
my $logfh;
my $snps_per_gene = {};
my $snp_nucleotides = {};
####################################################

my %options;
my $results = GetOptions (\%options,
						  "input_file|i=s",
						  "output_file|o=s",
						  "log|l=s",
						  "debug|d=s",
						  "help|h"
						 );

&check_options(\%options);

my $table = &parse_nucmer_snp_file( $options{'input_file'} );
$table->print_to_file( $options{'output_file'} );
print "$options{'output_file'}\n";

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

sub add_snp_counts {
  my ($table) = @_;
  foreach my $row ( $table->get_rows ) {

	## Add the query codon/aa
	unless( $row->gene_name eq 'intergenic' ) {
	  die("Could not find snp codon/aa info for ".$row->refpos) 
		unless( exists( $snp_nucleotides->{$row->refpos} ) );

	  my (@codons, @aas);
	  map {
		push(@codons, $snp_nucleotides->{$row->refpos}->{$_}->{'codon'} );
		push(@aas, $snp_nucleotides->{$row->refpos}->{$_}->{'aa'});
	  } keys %{$snp_nucleotides->{$row->refpos}};
	  $row->snp_codon( join("/", @codons) );
	  $row->snp_aa( join("/", @aas) );
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

  my $row = $table->get_row_by_position( $cols[10], $cols[0] );
  
  if ( !defined( $row ) ) {
	$row = new SNP::MergedTable::Row;

	# Check to see if there is an indel here.
	# merged table doesn't support indels.
	if ( $cols[1] eq '.' || $cols[2] eq '.' ) {
	  die("Input file should not have indels. Please use input file without indels");
	}

	$row->molecule( $cols[10] );
	$row->refpos( $cols[0] );
	$row->syn( $cols[15] );
	$row->refbase( $cols[1] );
	$row->gene_name( $cols[12] );
	$row->product( $cols[17] );
	$row->gene_start( $cols[13] );
	$row->gene_stop( $cols[14] );
	$row->pos_in_gene( $cols[16] );
	$row->ref_codon( $cols[19] );
	$row->ref_aa( $cols[20] );

	$row->query_base( $cols[11], $cols[2] );
	
	if( $row->gene_name eq 'intergenic' ) {
	  $row->gene_length( "NA" );
	  $row->snp_codon( "NA" );
	  $row->snp_aa( "NA" );
	} else {
	  ## If the query base is different from reference
	  if( $row->refbase ne $cols[2] ) {
		$snp_nucleotides->{$row->refpos} = {} unless( exists( $snp_nucleotides->{$row->refpos} ) );
		$snp_nucleotides->{$row->refpos}->{$cols[2]} = {'codon' => $cols[21],'aa' => $cols[22]};
	  }

	  my @lengths = ();
	  my @starts = split( "/", $cols[13] );
	  my @stops = split( "/", $cols[14] );
	  for ( my $i = 0; $i < @starts; $i++ ) {
		my $tmp = ($row->product eq 'NA') ? 'NA' : ($stops[$i] - $starts[$i]) + 1;
		push(@lengths, $tmp);
	  }
	  $row->gene_length( join("/", @lengths) );

	  ## Count snps_per_gene
	  # Sometimes we have a SNP located in multiple genes. Make sure to count it for both genes.
	  my @genes = split("/", $cols[12]);
	  foreach my $g ( @genes ) {
		$snps_per_gene->{$cols[10]} = {} unless( exists( $snps_per_gene->{$cols[10]} ) );
		unless ( exists( $snps_per_gene->{$cols[10]}->{$g} ) ) {
		  $snps_per_gene->{$cols[10]}->{$g} = 0;
		}
	  
		$snps_per_gene->{$cols[10]}->{$g}++;

	  }
	}
	
	$table->add_row( $row );

  } else {
	my $base = $row->query_base( $cols[11] );
	die("Already have query base for query $cols[11] [refpos: $cols[0]]") if( defined( $base ) );
	$row->query_base( $cols[11], $cols[2] );

	unless( $row->gene_name eq 'intergenic' ) {

	  ## If the query base is different from reference
	  if( $row->refbase ne $cols[2] ) {
		$snp_nucleotides->{$row->refpos} = {} unless( exists( $snp_nucleotides->{$row->refpos} ) );
		$snp_nucleotides->{$row->refpos}->{$cols[2]} = {'codon' => $cols[21],'aa' => $cols[22]};
	  }
	}
  }
  
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
