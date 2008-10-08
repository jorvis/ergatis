#!/usr/bin/perl

#./gff2fasta.pl [featuretype] [keysfile]
#./gff3_to_fasta.pl --type|t <featuretype> --input_file|i <file> --output_fasta|o <file>
# featuretype is stored in column 3 of gff3 
# For example, this is a line with featuretype cds
# nmpdr|158878.1.contig.NC_002758 NMPDR   cds ...
# only sequences of featuretype will be added to fasta file
use strict;
use warnings;

use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);

#use CGI qw/escape unescape/;
#use Data::Dumper;

my %opts = &parse_options();

my $type = $opts{type};

my $features = {};
my $fastadirective=0;
my $savefastaentry=0;
my $currfeature;

my $featcount = 0;
my $seqcount = 0;

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
    }
    else{
      #parse fasta
      if($line =~ /^>(\S+)/){
	$savefastaentry=0;
	if ( exists($features->{$1}) ) {
	  exists ($features->{$1}->{'save'}) || warn "save not defined for id ($1) on line ($line) dumped in ($opts{input_file}): ".Dumper($features->{$1});
	  if( $features->{$1}->{'save'} == 1) {
	    $savefastaentry=1;
	    ++$seqcount;
	    --$features->{$1}->{'save'};
	  }
	  elsif ($features->{$1}->{'save'} > 1){
	    die "Seen feature $1 ".$features->{$1}->{'save'}." > 1 times in $opts{input_file}";
	  }
	  else { # if save was decremented to 0, then we saw this already
	    die "Feature $1 associated with > 1 sequence";
	  }
	}
	else {
	  #print "No id ($1)\n";
	}
      }
      
      if($savefastaentry){
	print {$FOUT} $line;
      }
    }
  }
}

# if we never found a sequence for a feature, then save will still be 1
if ($featcount != $seqcount) {
  my $miss_seq_txt ='';
  foreach my $id (keys %{$features}) {
    $miss_seq_txt .= "missing sequence:\t$opts{input_file}\t$id\n" if ($features->{$id}->{save} == 1);
  }
  warn "Feature ($featcount) != Seqcount ($seqcount).  Missing features:\n".$miss_seq_txt;
}

# die if we never found any sequences
die "Missing all sequences in $opts{input_file}" if (-z $opts{output_fasta});



print "Feature count: $featcount\nSequence count: $seqcount\n";


#
# Subs
#

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
        ) || die "Unprocessable option";

    (defined($options{input_file})) || die "input_file is required parameter";
    # we may be passed a file name, but the .gz is what is actually there
    unless (-e $options{input_file}) {
    $options{input_file} .= ".gz";
    }
    (-r $options{input_file}) || die "input_file ($options{input_file}) not readable: $!";

    return %options;
}
