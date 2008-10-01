#!/usr/bin/perl

#./gff2fasta.pl [featuretype] [keysfile]
#./gff3_to_fasta.pl --type|t <featuretype> --input_gff3|i <file> --output_fasta|o <file>
# featuretype is stored in column 3 of gff3 
# For example, this is a line with featuretype cds
# nmpdr|158878.1.contig.NC_002758 NMPDR   cds ...
# only sequences of featuretype will be added to fasta file
use strict;
use warnings;

use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);

use CGI qw/escape unescape/;
use Data::Dumper;


my %opts = &parse_options();

my $type = $opts{type};

my $features = {};
my $fastadirective=0;
my $savefastaentry=0;
my $currfeature;
my @pepfastaoutbuffer;
#my @seqfastaoutbuffer;

my $featcount = 0;
my $seqcount = 0;

open (my $FIN, $opts{input_file}) || die "Unable to open input_file ($opts{input_file}): $!";
open (my $FOUT, ">$opts{output_fasta}") || die "Unable to open output_fasta ($opts{output_fasta}): $!";

#while(my $line=<STDIN>){
while(my $line=<$FIN>){
  $line = unescape($line);
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

	my @split_atts = split(/[;=]/,$elts[8]);
	
	die "No attributes on row" if (@split_atts == 0);
	die "Problem parsing attributes on row" if (@split_atts % 2 == 1);
	
	my(%attrs) =@split_atts;

	#exists($attrs{ID}) || die "Unable to parse ID from attributes";
	
	#die "line: $line" if ($attrs{ID} eq 'pathema|gb21.contig.76');

	#my(%dbxrefs) = split(/[,:]/,$attrs{'Dbxref'});
	$features->{$attrs{'ID'}}->{'save'}=1;

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
	  exists ($features->{$1}->{'save'}) || warn "save not defined for id ($1) on line ($line) dumped: ".Dumper($features->{$1});
	  if( $features->{$1}->{'save'} == 1) {
	    ++$features->{$1}->{'save'};
	    $savefastaentry=1;
	    ++$seqcount;
	  }
	  elsif ($features->{$1}->{'save'} > 1){
	    die "Seen feature $1 more than once";
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
    $miss_seq_txt .= "$id\n" if ($features->{$id}->{save} == 1);
  }
  die "Feature ($featcount) != Seqcount ($seqcount).  Missing features:\n".$miss_seq_txt;
}

print "Feature count: $featcount\nSequence count: $seqcount\n";


# if(defined $ARGV[1]){
#     open FILE,">>$ARGV[1]" or die;
#     foreach my $id (keys %$features){
# 	if($features->{$id}->{'type'} eq $ARGV[0]){
# 	    print FILE "$id $features->{$id}->{'center'},$features->{$id}->{'taxon'},$features->{$id}->{'genomic_source'},$features->{$id}->{'description'}\n";
# 	}
#     }
#     close FILE;
# }

# If(defined $ARGV[2]){
#     open FILE,">>$ARGV[2]" or die;
#     print FILE @seqfastaoutbuffer;
#     close FILE;
# }

# print @pepfastaoutbuffer;



#
# Subs
#

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
