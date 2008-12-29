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
      $line = clean_gff3_line($line);
      my (@elts) = split(/\t/,$line);
      $elts[2] = lc($elts[2]);

      # if it's a contig or supercontig (apidb), pull taxon info
      if($elts[2] eq 'contig' || $elts[2] eq 'source' || $elts[2] eq 'supercontig') {

	my(%attrs) = field_to_attributes($elts[8]);
	$attrs{ID} = exists($attrs{ID}) ? $attrs{ID} : $elts[0];

#	my(%dbxrefs) = split(/[,:]/,$attrs{'Dbxref'});

#	exists($dbxrefs{taxon}) || die "No taxon for contig $attrs{ID}";

	# try to parse out a taxon directly if one couldn't be pulled
	unless (  exists($attrs{taxon}) ) {
	  if ( $elts[8] =~ /taxon:(\d+)/i ) {
	    $attrs{taxon} = $1;
	  }
	}
	
	exists($attrs{taxon}) || die "No taxon for $attrs{ID} in $opts{input_file}".Dumper(%attrs);

	if ( (exists ($contig2taxon{$attrs{ID}})) && ($contig2taxon{$attrs{ID}} != $attrs{taxon})    ) {
	  # we will assume that this is for phage sequences and not really a problem
	  #      warn "Redef of taxon for $attrs{ID} ($opts{input_file}) ".$contig2taxon{$attrs{ID}}." != $dbxrefs{taxon}";
	  warn "Redef of taxon for $attrs{ID} ($opts{input_file}) ".$contig2taxon{$attrs{ID}}." != $attrs{taxon}";
	}
	else  {
	  $contig2taxon{$attrs{ID}} = $attrs{taxon};
	}

#	exists ($contig2taxon{$attrs{ID}}) && die "Redef of taxon for $attrs{ID} ".$contig2taxon{$attrs{ID}};
#	$contig2taxon{$attrs{ID}} = $dbxrefs{taxon};

      }

      if($elts[2] eq $type){

	my(%attrs) = field_to_attributes($elts[8]);
# 	my $id;
# 	if (defined( $attrs{ID} ) ) {
# 	  $id = $attrs{ID};
# 	}
# 	elsif (defined( $attrs{protein_id} )) {
# 	  $id = $attrs{protein_id};
# 	}
# 	#    unless (defined( $attrs{protein_id} )) {
# 	else {
# 	  # manually parse out protein_id
# 	  $elts[8] =~ q{protein_id=([^;]+)[;$]};
#           $attrs{protein_id}  = $1;
# 	  unless ( defined( $attrs{protein_id} )) {
# 	    warn "No protein_id for $attrs{ID} ($ifile)";
# 	    #next;
# 	  }
# 	}
	defined( $contig2taxon{$elts[0]} ) || die "No taxon for $elts[0] ($opts{input_file})";
#	defined($id) || die "no id";

	# make sure an ID was defined
	# can't go in field_to_attributes because it doesn't apply to contigs
	unless ( defined( $attrs{ID} ) ) {
	  # if there's no ID use locus_id or protein_id
	  if (defined( $attrs{locus_tag} )) {
	    $attrs{ID} = $attrs{locus_tag};
	  }
	  elsif (defined( $attrs{protein_id} )) {
	    $attrs{ID} = $attrs{protein_id};
	  }
	  else {
	    # manually parse out protein_id
	    if ($elts[8] =~ q{protein_id=([^;]+)[;$]} ) {
	       $attrs{ID} = $1;
	    }
	    else {
	       die "No protein_id for $attrs{ID} ($opts{input_file})";
	       next;
	    }
	  }
	}

        defined($attrs{ID}) || die "Missing 'ID' in attributes ( $opts{input_file} )";

	# unless a default source was passed in use column 2
	my $source = (exists $opts{source}) ? $opts{source} : $elts[1];
	my $contig = $elts[0];

	exists ($contig2taxon{$contig}) || die "Missing taxon for contig ($contig)";
	my $taxon = $contig2taxon{$contig};
       
	# check that a description was provided
	# if not use product
	unless ( defined($attrs{description})) {
	  if (defined( $attrs{product} )) {
	    $attrs{description} = $attrs{product};
	  }
	  else {
            # this really shouldn't be a die, as we know that not every cds will be annotated
	    warn "Missing 'description' in attributes ( $opts{input_file} )";
            $attrs{description} = '';
	  }
	}
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

  # remove trailing focus=
  $field =~ s/;focus=$//;

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
  
  my %attrs;
  
  while ((my $key = shift @split_atts) && (my $value = shift @split_atts)) {
    if ( ($key eq 'db_xref') || ($key eq 'Dbxref') ) {
      my @dbxrefs = split(/[:,]/,$value);
      while ((my $db = shift @dbxrefs) && (my $acc = shift @dbxrefs)) {
	$attrs{$db} = $acc;
      }
    }
    else {
      $attrs{$key} = $value;
    }
  }


  return %attrs;

}

sub clean_gff3_line {
  my $line = shift;

  # replace spaces
  $line =~ s/%20/ /g;

  # some files have them escaped, but just in the contig row, why??? :(
  $line =~ s/%7C/|/g;

  return $line;
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

