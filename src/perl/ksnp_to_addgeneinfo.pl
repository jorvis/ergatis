#!/usr/bin/perl

=head1 NAME

ksnp_merge_reference.pl 

B<--input,-in,-i>
    single fasta file or list of fasta files

=head1  DESCRIPTION

	Convert kSNP output to nucmer-show-snps output, allowing this information to be
	piped into the snp-add-gene-info component.

=head1  INPUT

	kSNP raw output with only one reference genome

=head1  OUTPUT

	The same as those output from nucmer-show-snps [see show-snps documentation for 
	description of columns]:
	p1, ref_base, query_base, p2, buff, dist, len_r, len_q, frm1, frm2, ref_contig, query_contig
    
=head1  CONTACT

    Kent Shefchek
    kshefchek@som.umaryland.edu

=cut

use strict;
use warnings;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Data::Dumper;

my $input;
my %options;
my $i = 0;
my $pkey = 0;
my %htable;
my @infile;
my $bool = 'true';

my $results = GetOptions(\%options,
			"input|in|i=s" => \$input
						 );

die "Must enter input file" unless defined $input;

open (my $fh, "<$input") || die "Cannot open input file";

#Split tab delimited kSNP output

@infile = <$fh>;

close $fh;

# determine orientation of file (order of query and reference varies based
# on alphabetical order

$bool= 'false' if ($infile[1] =~ /(.*)\t(.*)\t(.*)\t(\d+)\s(R|F)/);

foreach my $in (@infile){

$in =~ s/\n//g;

     if ($in !~ /^\s*$/){

     my ($kmerID, $kseq, $base, $pos, $fashead) = split(/\t/, $in);

     die "unexpected null value" if (($kmerID eq '') || ($kseq eq '') || ($base eq '') || ($pos eq '') || ($fashead eq ''));


          if ($kmerID != $i){

	  	if ($bool eq 'true'){
	  	&run_process_query_first(\%htable);
		}
		
		else {&run_process_ref_first(\%htable);}
	  
		
          #Empty table and reload first row

	   undef %htable;
	   $pkey = 0;

	   $htable{$pkey} = {
		        'kmerID' => $kmerID,
			'kseq' => $kseq,
			'base' => $base,
			'pos' => $pos,
			'fashead' => $fashead
		     };

   	        $pkey++; #increase primary key

	   $i++;

          }
	  else {

		#load hash table, hash contains a numeric keys pointing to hash references
 
     		$htable{$pkey} = {
			'kmerID' => $kmerID,
			'kseq' => $kseq,
			'base' => $base,
			'pos' => $pos,
			'fashead' => $fashead
		     };

    		$pkey++; #increase primary key
     		}
     }
}

#Run the subroutine one last time to get the last kmer group

if ($bool eq 'true'){
	  	&run_process_query_first(\%htable);
		}
		
		else {&run_process_ref_first(\%htable);}

##### Subroutines ##########

sub run_process_query_first{

my $hashref = shift;
my %table = %$hashref;
my ($qbase, $refbase, $refpos, $qfasta, $rfasta,$orient);

# Check to make sure that first row contains the query sequence

die "Unexpected value in query row" unless ($table{0}->{'pos'} eq 'x');

# Check to make sure that second row contains reference info 

die "Unexpected value in reference row\n".
    "Only use one reference genome for nucmer output" if ($table{1}->{'pos'} eq 'x') || ($table{1}->{'pos'} !~ /(\d+)\s(R|F)/);

$qfasta = $table{0}->{'fashead'};

my $len = keys %table;

for (my $i = 1; $i < $len; $i++){

     #reset query base if reverse complimented
 
     $qbase = $table{0}->{'base'}; 

     $rfasta = $table{$i}->{'fashead'};
     $refbase = $table{$i}->{'base'};
     $refpos = $table{$i}->{'pos'};

     ($refpos, $orient) = split (/\s/, $refpos);

     if ($orient eq 'R'){

	$refbase = &reverse_complement($refbase);
        $qbase = &reverse_complement($qbase);

     }

     # print p1, ref_base, query_base, p2, buff, dist, len_r, len_q, frm1, frm2, ref_contig, query_contig

     print STDOUT "$refpos\t$refbase\t$qbase\t0\t0\t0\t0\t0\t1\t0\t$rfasta\t$qfasta\n";


}

}

sub run_process_ref_first{

my $hashref = shift;
my %table = %$hashref;
my ($qbase, $refbase, $refpos, $qfasta, $rfasta,$orient);

# Check to make sure that first row contains the query sequence

my $len = keys %table;

my $lastkey = $len - 1;

die "Unexpected value in query row" unless ($table{$lastkey}->{'pos'} eq 'x');

# Check to make sure that second row contains reference info 

die "Unexpected value in reference row\n".
    "Only use one reference genome for nucmer output" if ($table{0}->{'pos'} eq 'x') || ($table{0}->{'pos'} !~ /(\d+)\s(R|F)/);

$qfasta = $table{$lastkey}->{'fashead'};
$qbase = $table{$lastkey}->{'base'};
my $cubase= $qbase;

delete $table{$lastkey};

$len = keys %table;

for (my $i = 0; $i < $len; $i++){

     $qbase = $cubase;

     $rfasta = $table{$i}->{'fashead'};
     $refbase = $table{$i}->{'base'};
     $refpos = $table{$i}->{'pos'};

     ($refpos, $orient) = split (/\s/, $refpos);

     if ($orient eq 'R'){

	$refbase = &reverse_complement($refbase);
        $qbase = &reverse_complement($qbase);

     }

     # print p1, ref_base, query_base, p2, buff, dist, len_r, len_q, frm1, frm2, ref_contig, query_contig

     print STDOUT "$refpos\t$refbase\t$qbase\t0\t0\t0\t0\t0\t1\t0\t$rfasta\t$qfasta\n";


}

}


sub reverse_complement {
	my ($seq) = @_;
	my $rev = reverse $seq;
	$rev =~ tr/ATCGRYWSKMVBHDatcgrywskmvbhd/TAGCYRWSMKBVDHtagcyrwsmkbvdh/;
	return $rev;
}


















