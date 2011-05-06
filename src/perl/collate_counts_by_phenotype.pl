#!/usr/bin/perl

=head1 NAME

collate_couns_by_phenotype.pl - given a file with a list of read count files
and a file with how filename prefixes in those files match to a replicate and phenotyp group,
this script will output a tab-delimeted file like this:
<replicate id> <phenotype> <read count file>

The file names in the read count files must be unique so that they can be matched to replicates 
using the information in the sample matching file.



=head1 SYNOPSIS

  USAGE: collate_couns_by_phenotype.pl <read count list file> <sample matching file> <output file>

=cut

use strict;

my $count_paths = $ARGV[0];
my $sample_matching = $ARGV[1];
my $output = $ARGV[2];

my %replicates;

open COUNTS, "$count_paths";
my @COUNT_FILES = <COUNTS>;
close COUNTS;

open REPLICATES, "$sample_matching";

while (<REPLICATES>){
    chomp $_;
    my ($replicate, $phenotype, $pattern) = split("\t", $_);
    
    foreach(@COUNT_FILES){

	chomp $_;

	my $filename;

	if ($_ =~ /.*\/(.+)$/){
	    $filename = $1;
	}
	else{
	    die "\nError! unable to get basename from $_!";	    
	}
	
	if($filename =~ /$pattern/){
	    $replicates{$replicate} = { 'phenotype' => $phenotype, 'pattern' => $pattern, 'count_file' => $_ };
	    last;
	}
    }
}

close REPLICATES;

open OUT, ">$output";

foreach (sort keys %replicates){
    print OUT join("\t", $_, $replicates{$_}{'phenotype'}, $replicates{$_}{'count_file'})."\n";
}

close OUT;







	
	

	
