#!/usr/bin/perl

#./gff3_to_annotab.pl --type|t <featuretype> --input_file|i <file> --output_annotab|o <file> [--source|s <source>]
# featuretype is stored in column 3 of gff3 
# For example, this is a line with featuretype cds
# nmpdr|158878.1.contig.NC_002758 NMPDR   cds ...
# output is tab-delimited file with columns
#cds_id    source    taxon_id    contig_id    product_name    [GO assignments]    [EC numbers]
# ie
# pathema|ntbc02.1.cds.ORF00001   Pathema 288681  pathema|ntbc02.contig.1 chromosomal replication initiator protein DnaA
use strict;
use warnings;

use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);

#use CGI qw/escape unescape/;
#use Data::Dumper;

my %opts = &parse_options();

my $type = $opts{type};

my $fastadirective=0;

my %contig2taxon;

my $featcount = 0;
my $seqcount = 0;

open (my $FIN, $opts{input_file}) || die "Unable to open input_file ($opts{input_file}): $!";
open (my $FOUT, ">$opts{output_annotab}") || die "Unable to open output_annotab ($opts{output_annotab}): $!";

while(my $line=<$FIN>){
  #$line = unescape($line);
  if($line =~ /^\#\#FASTA/){
    $fastadirective=1;
  }
  if($line !~ /^\#/){
    if(!$fastadirective){
      #parse tab delimeted
      chomp $line;
      my (@elts) = split(/\t/,$line);
      $elts[2] = lc($elts[2]);

      # if it's a contig or supercontig (apidb), pull taxon info
      if($elts[2] eq 'contig' || $elts[2] eq 'supercontig') {

	my(%attrs) = field_to_attributes($elts[8]);

	my(%dbxrefs) = split(/[,:]/,$attrs{'Dbxref'});

	exists($dbxrefs{taxon}) || die "No taxon for contig $attrs{ID}";

	exists ($contig2taxon{$attrs{ID}}) && die "Redef of taxon for $attrs{ID} ".$contig2taxon{$attrs{ID}};
	$contig2taxon{$attrs{ID}} = $dbxrefs{taxon};

      }

      if($elts[2] eq $type){

	my(%attrs) = field_to_attributes($elts[8]);

	# unless a default source was passed in use column 2
	my $source = (exists $opts{source}) ? $opts{source} : $elts[1];
	my $contig = $elts[0];

	exists ($contig2taxon{$contig}) || die "Missing taxon for contig ($contig)";
	my $taxon = $contig2taxon{$contig};
       
	defined($attrs{ID}) || die "Missing 'ID' in attributes ( $opts{input_file} )";
	defined($attrs{description}) || die "Missing 'description' in attributes ( $opts{input_file} )";

	print {$FOUT} join("\t", ($attrs{ID}, $source, $contig2taxon{$contig}, $contig, $attrs{description}));

	# optional fields

	if ($opts{parse_GO}) {
	  print {$FOUT} "\t".parse_GO(%attrs); # null or comma-delimited GO ids
	}

	if ($opts{parse_EC}) {
	  print {$FOUT} "\t".parse_EC(%attrs); # null or single EC id
	}

	if ($opts{parse_gene_symbol}) {
	  print {$FOUT} "\t".parse_gene_symbol(%attrs); # null or value of gene_symbol
	}

	print {$FOUT} "\n";
	++$featcount;
      }
    }
  }
}

#die "Feature ($featcount) != Seqcount ($seqcount)" if ($featcount != $seqcount);
print "Feature count: $featcount\n";

#
# Subs
#

# Ways EC numbers might be stored:
# pathema: Dbxref=TIGR_CMR:BAKB_0005,EC:2.7.8.-;
sub parse_EC {
  my %attrs = @_;

  if (defined ($attrs{Dbxref})) {
    # assume comma delimited 
    # it looks like EC only matches once, but this is a bug if it matches > 1
    foreach my $full_dbxref (split(/[,]/,$attrs{'Dbxref'})) {
      if ($full_dbxref =~ /^EC:/) {
	return $full_dbxref;
      }
    }
  }
  return 'null';
}

# Ways GO ids might be stored:
# pathema:  Ontology_term=GO:0004239,GO:0006464
sub parse_GO {
  my %attrs = @_;

  if (defined ($attrs{Ontology_term})) {
    return $attrs{Ontology_term};
  } else {
    return 'null';
  }
}

# Ways gene_symbol might be stored:
# pathema:  gene_symbol=rpmI;
sub parse_gene_symbol {
  my %attrs = @_;

  if (defined ($attrs{gene_symbol})) {
    return $attrs{gene_symbol};
  } else {
    return 'null';
  }
}

#shared w/ gff3_to_fasta.pl
sub field_to_attributes {
  my $field = shift;

  #remove empty keypairs
  $field =~ s/;;+/;/g;

  #remove trailing keypairs
  $field =~ s/;\s*$//g;

  my @split_atts = split(/[;=]/,$field);

  # remove leading spaces
  foreach (@split_atts) {
    $_ =~ s/^\s+//;
  }

  # because above dies on ;; in attribute field...
#  foreach my $keypair ( split( /[;]/, $field) ) {
#    print "$keypair\n";
#  }
  
  die "No attributes on row" if (@split_atts == 0);
  die "Odd number of keys parsed from attribute field ($field) in (  $opts{input_file} )" if (@split_atts % 2 == 1);
  
  my(%attrs) =@split_atts;

  return %attrs;

}


sub parse_options {
    my %options = ();
    GetOptions( \%options,
        'input_file|i=s',
        'output_annotab|o=s',
        'type|t=s',
        'source|s=s',
	'parse_GO:i', # optional default 0
	'parse_EC:i',
	'parse_gene_symbol:i',	       
        ) || die "Unprocessable option";

    (defined($options{input_file})) || die "input_file is required parameter";
    # we may be passed a file name, but the .gz is what is actually there
#    unless (-e $options{input_file}) {
#    $options{input_file} .= ".gz";
#    }
    (-r $options{input_file}) || die "input_file ($options{input_file}) not readable: $!";

    return %options;
}
