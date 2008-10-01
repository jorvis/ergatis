#!/usr/bin/perl

#./gff3_to_annotab.pl --type|t <featuretype> --input_file|i <file> --output_annotab|o <file> [--source|s <source>]
# featuretype is stored in column 3 of gff3 
# For example, this is a line with featuretype cds
# nmpdr|158878.1.contig.NC_002758 NMPDR   cds ...
use strict;
use warnings;

use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);

use CGI qw/escape unescape/;
#use Data::Dumper;

my %opts = &parse_options();

my $type = $opts{type};

#my $features = {};
my $fastadirective=0;
#my $savefastaentry=0;
#my $currfeature;
#my @pepfastaoutbuffer;
#my @seqfastaoutbuffer;

my %contig2taxon;

my $featcount = 0;
my $seqcount = 0;

open (my $FIN, $opts{input_file}) || die "Unable to open input_file ($opts{input_file}): $!";
open (my $FOUT, ">$opts{output_annotab}") || die "Unable to open output_annotab ($opts{output_annotab}): $!";

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

      # if it's a contig, pull taxon info
      if($elts[2] eq 'contig') {
	my @split_atts = split(/[;=]/,$elts[8]);
	
	die "No attributes on row" if (@split_atts == 0);
	die "Problem parsing attributes on row" if (@split_atts % 2 == 1);
	
	my(%attrs) =@split_atts;

	my(%dbxrefs) = split(/[,:]/,$attrs{'Dbxref'});

	exists($dbxrefs{taxon}) || die "No taxon for contig $attrs{ID}";

	exists ($contig2taxon{$attrs{ID}}) && die "Redef of taxon for $attrs{ID} ".$contig2taxon{$attrs{ID}};
	$contig2taxon{$attrs{ID}} = $dbxrefs{taxon};

      }

      if($elts[2] eq $type){

	my @split_atts = split(/[;=]/,$elts[8]);
	
	die "No attributes on row" if (@split_atts == 0);
	die "Problem parsing attributes on row" if (@split_atts % 2 == 1);
	
	my(%attrs) =@split_atts;

	# unless a default source was passed in use column 2
	my $source = (exists $opts{source}) ? $opts{source} : $elts[1];
	my $contig = $elts[0];

	exists ($contig2taxon{$contig}) || die "Missing taxon for contig ($contig)";
	my $taxon = $contig2taxon{$contig};
       
	defined($attrs{ID}) || die "Missing 'ID' in attributes";
	defined($attrs{description}) || die "Missing 'description' in attributes";

	print {$FOUT} join("\t", ($source, $contig, $contig2taxon{$contig}, $attrs{ID}, $attrs{description}))."\n";

	++$featcount;
      }
    }
#     else{
#       #parse fasta
#       if($line =~ /^>(\S+)/){
# 	$savefastaentry=0;
# 	if ( exists($features->{$1}) ) {
# 	  exists ($features->{$1}->{'save'}) || warn "save not defined for id ($1) on line ($line) dumped: ".Dumper($features->{$1});
# 	  if( $features->{$1}->{'save'} == 1) {
# 	    ++$features->{$1}->{'save'};
# 	    $savefastaentry=1;
# 	    ++$seqcount;
# 	  }
# 	  elsif ($features->{$1}->{'save'} > 1){
# 	    die "Seen feature $1 more than once";
# 	  }
# 	}
# 	else {
# 	  #print "No id ($1)\n";
# 	}
#       }
      
#       if($savefastaentry){
# 	print {$FOUT} $line;
#       }
#     }
  }
}

#die "Feature ($featcount) != Seqcount ($seqcount)" if ($featcount != $seqcount);
print "Feature count: $featcount\n";


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
        'output_annotab|o=s',
        'type|t=s',
        'source|s=s',
        ) || die "Unprocessable option";

    (defined($options{input_file})) || die "input_file is required parameter";
    # we may be passed a file name, but the .gz is what is actually there
    unless (-e $options{input_file}) {
    $options{input_file} .= ".gz";
    }
    (-r $options{input_file}) || die "input_file ($options{input_file}) not readable: $!";

    return %options;
}
