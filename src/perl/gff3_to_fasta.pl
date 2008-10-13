#!/usr/bin/perl

#./gff2fasta.pl [featuretype] [keysfile]
#./gff3_to_fasta.pl --type|t <featuretype> --input_file|i <file> --output_fasta|o <file> [--sequence_type|s <featuretype]
# featuretype is stored in column 3 of gff3 
# For example, this is a line with featuretype cds
# nmpdr|158878.1.contig.NC_002758 NMPDR   cds ...
# only sequences of type featuretype will be added to fasta file
# the optional parameter sequence_type specifies the type of the child features assiated
# with the Parent featuretype that should be used for the Fasta sequence
# This is needed to associate cds sequences with their parent mrnas

use strict;
use warnings;

use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);

#use CGI qw/escape unescape/;
use Data::Dumper;

my %opts = &parse_options();

my $type = $opts{type};
my $sequence_type = $opts{sequence_type};

my $features = {};
my $fastadirective=0;
my $savefastaentry=0;
my $currfeature;

my $featcount = 0;
my $seqcount = 0;

# map a parent id to its children
# used for linking sequence_type features to the main featuretype
my %parent2id;
my %id2parent;

open (my $FIN, $opts{input_file}) || die "Unable to open input_file ($opts{input_file}): $!";
open (my $FOUT, ">$opts{output_fasta}") || die "Unable to open output_fasta ($opts{output_fasta}): $!";

while(my $line=<$FIN>){
#  $line = unescape($line);
  if($line =~ /^\#\#FASTA/){
    $fastadirective=1;
  }
  if($line !~ /^\#/){
    if(!$fastadirective){
      #parse tab delimeted
      chomp $line;
      my (@elts) = split(/\t/,$line);
      $elts[2] = lc($elts[2]);

      if($elts[2] eq $type){

	my(%attrs) = field_to_attributes($elts[8]);

	#exists($attrs{ID}) || die "Unable to parse ID from attributes";
	
	#my(%dbxrefs) = split(/[,:]/,$attrs{'Dbxref'});
	++$features->{$attrs{'ID'}}->{'save'};

	++$featcount;

	# only care about save, note that taxon causes problems
#	$features->{$attrs{'ID'}}->{'description'} = $attrs{'description'};
#	$features->{$attrs{'ID'}}->{'center'} = $elts[1];
#	$features->{$attrs{'ID'}}->{'genomic_source'} = $elts[0];
#	$features->{$attrs{'ID'}}->{'taxon'} = $features->{$elts[0]}->{'taxon'};
#	$features->{$attrs{'ID'}}->{'type'} = $elts[2];

      }
      elsif ($elts[2] eq $sequence_type){ # hopefully we never get a sequence type of 0
	my(%attrs) = field_to_attributes($elts[8]);

	foreach ('ID', 'Parent') {
	  exists($attrs{$_}) || die "Unable to parse $_ from attributes";
	}

	++$parent2id{$attrs{Parent}}->{$attrs{ID}}; # increase count of this map
	++$id2parent{$attrs{ID}}->{$attrs{Parent}}; # increase count of this map

      }
    }
    else{
      #parse fasta
      if($line =~ /^>(\S+)/){
	$savefastaentry=0;
	my $id = $1;
	# are we outputting child sequences instead of the parent?
	if ( $sequence_type ) {
	  next unless (defined $id2parent{$id});
	  if ( scalar(keys %{$id2parent{$id}}) > 1 ) {
	    die "Child id $id mapped to > 1 Parents:\t". join("\t",(keys %{$id2parent{$id}}));
	  }

	  my $parent = (keys %{$id2parent{$id}})[0]; # grab the parent id
	  
	  if ( scalar(keys %{$parent2id{$parent}}) > 1 ) {
	    die "Parent id $parent mapped to > 1 Childen:\t". join("\t",(keys %{$parent2id{$parent}}));
	  }

	  # if the parent was of featuretype update seen counts and set savefastaentry=1
	  if ( exists($features->{$parent}) ) {
	    track_features_printed($parent, $line);
	    $line =~ s/^>$id/>$parent\t$id/; # update fasta header
	  }
	}
	# if we don't care about parent type, check if this is of featuretype
	elsif ( exists($features->{$id}) ) {
	  track_features_printed($id, $line);
	}
	else {
	  #print "No id ($id)\n";
	}
      }
      
      if($savefastaentry){
	print {$FOUT} $line;
      }
    }
  }
}

close($FOUT);

# if we never found a sequence for a feature, then save will still be 1
if ($featcount != $seqcount) {
  my $miss_seq_txt ='';
  foreach my $id (keys %{$features}) {
#    $miss_seq_txt .= "missing sequence:\t$opts{input_file}\t$id\n" if ($features->{$id}->{save} == 1);
    $miss_seq_txt .= "\t$id" if ($features->{$id}->{save} == 1);
  }
  warn "Feature / sequence count mismatch ($featcount != $seqcount). Missing features: ".$miss_seq_txt,"\n";
}

# die if we never found any sequences
die "Empty fasta output from $opts{input_file}" if (-z $opts{output_fasta});



print "Feature count: $featcount\nSequence count: $seqcount\n";


#
# Subs
#

# check if a feature is one that we want to return
# and if its been returned the correct number of times
sub track_features_printed {
  my $id = shift;
  my $line = shift;
	  exists ($features->{$id}->{'save'}) || warn "save not defined for id ($id) on line ($line) dumped in ($opts{input_file}): ".Dumper($features->{$id});
	  if( $features->{$id}->{'save'} == 1) {
	    $savefastaentry=1;
	    ++$seqcount;
	    --$features->{$id}->{'save'};
	  }
	  elsif ($features->{$id}->{'save'} > 1){
	    die "Seen feature $id ".$features->{$id}->{'save'}." > 1 times in $opts{input_file}";
	  }
	  else { # if save was decremented to 0, then we saw this already
	    die "Feature $id associated with > 1 sequence";
	  }

}

#shared w/ gff3_to_annotab.pl
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
        'output_fasta|o=s',
        'type|t=s',
        'sequence_type|s=s',
        ) || die "Unprocessable option";

    (defined($options{input_file})) || die "input_file is required parameter";
    # we may be passed a file name, but the .gz is what is actually there
    unless (-e $options{input_file}) {
    $options{input_file} .= ".gz";
    }
    (-r $options{input_file}) || die "input_file ($options{input_file}) not readable: $!";

    (defined($options{sequence_type})) || ($options{sequence_type} = 0); # optional, if not provided skip over dealing with it

    return %options;
}
