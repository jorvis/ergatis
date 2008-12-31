#!/usr/bin/perl

#./gff3_to_annotab.pl --type|t <featuretype> --input_file|i <file> --output_annotab|o <file> [--source|s <source>] [--id_attr <attr_field>]
# seqID will be taken from id_attr if one is specified.  If id_attr is specified and a feature of <featuretype> lacks it, the script will die.
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
use Data::Dumper;

my %opts = &parse_options();

my $type = $opts{type};

my $fastadirective=0;

my %contig2taxon;

my $featcount = 0;

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
	$attrs{seqID} = exists($attrs{ID}) ? $attrs{ID} : $elts[0];

#	my(%dbxrefs) = split(/[,:]/,$attrs{'Dbxref'});

#	exists($dbxrefs{taxon}) || die "No taxon for contig $attrs{ID}";

	# try to parse out a taxon directly if one couldn't be pulled
	unless (  exists($attrs{taxon}) ) {
	  if ( $elts[8] =~ /taxon:(\d+)/i ) {
	    $attrs{taxon} = $1;
	  }
	}
	
	exists($attrs{taxon}) || die "No taxon for $attrs{seqID} in $opts{input_file}".Dumper(%attrs);

	if ( (exists ($contig2taxon{$attrs{seqID}})) && ($contig2taxon{$attrs{seqID}} != $attrs{taxon})    ) {
	  # we will assume that this is for phage sequences and not really a problem
	  #      warn "Redef of taxon for $attrs{ID} ($opts{input_file}) ".$contig2taxon{$attrs{ID}}." != $dbxrefs{taxon}";
	  warn "Redef of taxon for $attrs{seqID} ($opts{input_file}) ".$contig2taxon{$attrs{seqID}}." != $attrs{taxon}";
	}
	else  {
	  $contig2taxon{$attrs{seqID}} = $attrs{taxon};
	}

#	exists ($contig2taxon{$attrs{seqID}}) && die "Redef of taxon for $attrs{seqID} ".$contig2taxon{$attrs{seqID}};
#	$contig2taxon{$attrs{seqID}} = $dbxrefs{taxon};

      }

      if($elts[2] eq $type){

	my(%attrs) = field_to_attributes($elts[8]);

	defined( $contig2taxon{$elts[0]} ) || die "No taxon for $elts[0] ($opts{input_file})";

	# retrieve seqID from specific attribute field(s) if one was specified
	if ( defined($opts{id_attr}) && !($opts{id_attr} eq '')) {
	  my $found_id = 0;
	  my @id_attr = split(",", $opts{id_attr});
	  foreach my $id_attr (@id_attr) {
	    if ( defined($attrs{$id_attr}) ) {
	      $attrs{seqID} = $attrs{$id_attr};
	      ++$found_id;
	      last;
	    }
	  }
	  unless ($found_id) {
	    die "No specified id_attr $opts{id_attr} for row in file $opts{input_file}\n".Dumper(%attrs);
	  }
	}
	# the below is all largely redundant after comma-delimited method above
	# assume that if an id_attr was specified, then for consistency's sake all seqIDs must come from that
	else {
	  # make sure an seqID was defined
	  # can't go in field_to_attributes because it doesn't apply to contigs
	  my $found_id = 0;
	  foreach ('ID', 'locus_tag', 'protein_id') {
	    if ( defined($attrs{$_}) ) {
	      $attrs{seqID} = $attrs{$_};
	      ++$found_id;
	      last;
	    }
	  }
	  unless ($found_id) {
	    die "No specified id_attr $opts{id_attr} for row in file $opts{input_file}\n".Dumper(%attrs);
	    # manually parse out protein_id
	    if ($elts[8] =~ q{protein_id=([^;]+)[;$]} ) {
	       $attrs{seqID} = $1;
	    }
	    else {
	       die "Unable to populate seqID for line ($line) in $opts{input_file}";
	       next;
	    }
	  }
# 	unless ( defined( $attrs{ID} ) ) {
# 	  # if there's no ID use locus_id or protein_id
# 	  if (defined( $attrs{locus_tag} )) {
# 	    $attrs{ID} = $attrs{locus_tag};
# 	  }
# 	  elsif (defined( $attrs{protein_id} )) {
# 	    $attrs{ID} = $attrs{protein_id};
# 	  }
# 	  else {
# 	    # manually parse out protein_id
# 	    if ($elts[8] =~ q{protein_id=([^;]+)[;$]} ) {
# 	       $attrs{ID} = $1;
# 	    }
# 	    else {
# 	       die "No id for ($line) in $opts{input_file}";
# 	       next;
# 	    }
# 	  }
# 	}
        }

        defined($attrs{seqID}) || die "Unable to parse seqID from attributes ( $opts{input_file} )";

	# unless a default source was passed in use column 2
	my $source = (exists $opts{source}) ? $opts{source} : $elts[1];
	my $contig = $elts[0];

	exists ($contig2taxon{$contig}) || die "Missing taxon for contig ($contig)";
	my $taxon = $contig2taxon{$contig};
       
	# check that a description was provided
	# if not use product
	if ( !( defined($attrs{description}) ) && (defined($attrs{product})) ) {
	  $attrs{description} = $attrs{product};
	}

	if ( defined($attrs{description}) ) {
	  # remove all unicode characters
	  $attrs{description} =~ s/%([a-fA-F0-9]{2,2})/chr(hex($1))/eg;
	}
	else {
	  # this really shouldn't be a die, as we know that not every cds will be annotated
	  warn "Missing 'description' in attributes ( $opts{input_file} )";
	  $attrs{description} = '';
	}
	print {$FOUT} join("\t", ($attrs{seqID}, $source, $contig2taxon{$contig}, $contig, $attrs{description}));

	# optional fields
	if ($opts{parse_GO}) {
	  print {$FOUT} "\t".parse_GO($elts[8], %attrs); # null or comma-delimited GO ids
	}

	if ($opts{parse_EC}) {
	  print {$FOUT} "\t".parse_EC($elts[8], %attrs); # null or single EC id
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

print "Feature count: $featcount\n";

#
# Subs
#

# Ways EC numbers might be stored:
# pathema: Dbxref=TIGR_CMR:BAKB_0005,EC:2.7.8.-;
sub parse_EC {
  my $field = shift;
  my %attrs = @_;

  if ( defined($attrs{EC_number}) ) {
    return "EC:".$attrs{EC_number};
  }
  elsif ( defined($attrs{EC}) ) {
    return "EC:".$attrs{EC};
  }
  elsif ( $attrs{description} =~ /EC:[\.\-\d]+/ ) { #regex works w/o ()'s
    return join(",", sort ( $attrs{description} =~ /EC:[\.\-\d]+/g ) );
  }
  elsif ( $field =~ /EC:[\.\-\d]+/ ) {
    return join(",", sort ( $field =~ /EC:[\.\-\d]+/g ) );	      
  }

  return 'null';
}

# Ways GO ids might be stored:
# pathema:  Ontology_term=GO:0004239,GO:0006464
sub parse_GO {
  my $field = shift;
  my %attrs = @_;

  if (defined ($attrs{Ontology_term})) {
    return $attrs{Ontology_term};
  } 
  # GO_component:%200005840 or GO_process:%200009116 
  # but it looks like the characters have been decoded
  elsif ( defined($attrs{note}) && ( $attrs{note} =~ m/GO_[a-zA-Z]+:\s*(\d+)/ ) ) {
#  elsif ( defined($attrs{note}) && ( $attrs{note} =~ m/(GO_[a-zA-z]+:....)/ ) ) {
    my $go = "GO:$1";
#    $go  =~ s/%([a-fA-F0-9]{2,2})/chr(hex($1))/eg; # remove encodings
    return $go;
  }
  elsif (  defined($attrs{GO}) ) {
    return "GO:".$attrs{GO};
  }
  elsif ( $attrs{description} =~ /GO:[\d]+/ ) {
    return join(",", sort ( $attrs{description} =~ /GO:[\d]+/g ) );
  }
  elsif ( $field =~ /GO:[\d]+/ ) {
    return join(",", sort ( $field =~ /GO:[\d]+/g ) );	      
  }
  else {
    return 'null';
  }
}

# Ways gene_symbol might be stored:
# pathema:  gene_symbol=rpmI;
sub parse_gene_symbol {
  my %attrs = @_;

  if (defined ($attrs{gene_symbol})) {
    return $attrs{gene_symbol};
  } 
  # : separated portion of the ID, e.g.,
  # ID=CP000082.1:dnaN:unknown_transcript_1;
  elsif ( defined($attrs{ID}) && (my @id_parts = split(":", $attrs{ID})) > 1 ) {
#    die Dumper(%attrs);    
    return $id_parts[1];
  }
  else {
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

  # to handle empty valued keys (e.g., key=;)
  %attrs=@split_atts;
#  print "split_atts: \n".Dumper(@split_atts);
#  print "attrs:\n".Dumper(%attrs);

  foreach my $key ('db_xref', 'Dbxref') {
    if ( exists($attrs{$key}) ) {
      my @dbxrefs = split(/[:,]/,$attrs{$key});
      while ((my $db = shift @dbxrefs) && (my $acc = shift @dbxrefs)) {
	if ( defined($attrs{$db}) && ($attrs{$db} ne $acc) ) {
	  die "db value conflict for $db: $acc ne $attrs{$db}";
	}
 	$attrs{$db} = $acc;
      }
    }
  }

  return %attrs;

}

sub clean_gff3_line {
  my $line = shift;

  # replace spaces
  $line =~ s/%20/ /g;

  # replace (
  $line =~ s/%28/\(/g;

  # replace )
  $line =~ s/%29/\)/g;

  # some files have them escaped, but just in the contig row, why??? :(
  $line =~ s/%7C/|/g;

  # can't do this here because it'll decode ; (%3B) and screw up parsing
#  $line =~ s/%([a-fA-F0-9]{2,2})/chr(hex($1))/eg;

  return $line;
}

sub parse_options {
    my %options = ();
    GetOptions( \%options,
        'input_file|i=s',
        'output_annotab|o=s',
        'type|t=s',
        'source|s=s',
        'id_attr:s',
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

