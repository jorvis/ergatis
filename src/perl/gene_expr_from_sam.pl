#!/usr/bin/perl

=head1 NAME

gene_expr_from_sam.pl - parses a SAM alignment file and calculates expression data from it

=head1 SYNOPSIS

 USAGE: gene_expr_from_sam.pl sam_file

=head1  CONTACT

    Umar Farooq
    ufarooq@som.umaryland.edu

=cut

use strict;
use warnings;
use Pod::Usage;
use Data::Dumper;

open (SAM, $ARGV[0]) || die "can't open SAM file: $ARGV[0]!";

my %genes;
my $reads = 0;

while (my $line = <SAM>){
    chomp $line;
    if($line =~ /^\@SQ/){
	if($line=~ /SN:(\w+)\s+LN:(\d+)/){
	    my $len = $2 ? $2 : 0;
	    $genes{ $1 } = {'length' => $len, 'aligns' => 0, 'rpkm' => 0};
	    print "\n$1\t$genes{$1}{'length'}";
	}
	
    }   
    else{
	my @cols = split("\t", $line);
	if($cols[1] ne '77' && $cols[1] ne '141'){
	    if(@cols>8 && $cols[6] eq '='){
		if($cols[1] eq '163' ||	$cols[1] eq '83' || $cols[1] eq '99' ||	$cols[1] eq '147'){
		    #add half for each read, because we assume paired reads
		    $genes {$cols[2] } {'aligns'} += 0.5 ;
		    $reads += 0.5;
		}
		else {
		    die "some kinda problem with paired read $line";
		}
	    }
	}
    }
}

close SAM;

open OUT, ">$ARGV[0].expr";

print "\nTOTAL MAPPED READS: $reads";

foreach my $gene (sort(keys (%genes))){
    
    if(($genes{$gene}{'length'}>0) && ($reads>0)){
	$genes{$gene}{'rpkm'} = ($genes{$gene}{'aligns'}/($genes{$gene}{'length'}/1000))/($reads/1000000);
    }
    print OUT "\n$gene\t$genes{$gene}{'length'}\t$genes{$gene}{'aligns'}\t$genes{$gene}{'rpkm'}";
}

close OUT;

