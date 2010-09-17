#!/usr/bin/perl

=head1 NAME

gene_coverage_stats.pl - Parses data from the coverage output from MosaikCoverage,
                         and the assembler output from MosaikAssembler and prints a tab-delimited output file
                         that includes coverage information for each gene, as well as
                         Mapped Reads per Kb of Ref Region per million mapped reads (RPKM) for each locus

=head1 SYNOPSIS

 USAGE: gene_coverage_stats.pl
       --coverage_file=/path/to/coverage/file/output/from/mosaikcoverage
       --assembler_file=/path/to/assembler/file/output/from/mosaikassemble
       --output_file=/path/to/output/file
     [ --help
     ]

=head1 OPTIONS

B<--output_file>
    REQUIRED. Output file where output using the cutoffs is placed in tab-delimited format.  

B<--coverage_file>
    REQUIRED. Coverage file output from using MosaikBuild

B<--assembler_file>
    REQUIRED. Assembler file output from MosaikAssembler

B<--help,-h>
    Print this message

=head1  CONTACT

    Umar Farooq
    ufarooq@som.umaryland.edu

=cut

use strict;
use warnings;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Pod::Usage;
use Data::Dumper;

my %options = ();
my $results = GetOptions (\%options, 
						  'coverage_file=s',
			                          'assembler_file=s',
			                          'sort_file=s',
						  'output_file=s',
                          'help|h') || pod2usage();

# display documentation
if( $options{'help'} ){
    pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDERR} );
}


# make sure parameters are correct
&check_parameters(\%options);


my %cov;
my $total_mapped_reads = 0;

open IN, "$options{coverage_file}";
while( my $line = <IN> ) {
	chomp $line;

	next unless( $line =~ /^\*/ );

	my ($locus, $percentage, $mean, $gene_length);
	my (@items) = split(" ", $line);

	$locus = $items[4];

	if($items[5] =~ /\((\d+)/)
	    {
		$gene_length = $1;
	    }
	if($items[7]==0) {
	    $percentage = 0;
	    $mean = 0;
	}
	else {
	    $items[9] =~ /\((.+)/;
	    $percentage = $1;
	    $items[12] =~ /(.+)x/;
	    $mean = $1;
	}
	$cov{ $locus } = { 'gene_length' => $gene_length, 'coverage_length' => $items[7], 'percentage' => $percentage, 'mean_coverage' => $mean, 'RPKM' => 0, 'mapped_reads' => 0 };
}

close IN;

#read in the number of alignments for each gene from the assembler file
my $afile = $options{assembler_file};
open (AFILE, "$afile");
my $line = <AFILE>;
chomp $line;

do {
    $line = <AFILE>;
    chomp $line;
} until ($line eq "alignment count   reference sequence");
$line = <AFILE>;
chomp $line;
$line = <AFILE>;
chomp $line;

while($line =~ /^\s+(\d+)\s+(.+)/){
    my $locus = $2;

    #this will catch the case of MosaikAssemlber mentioning a gene not in the MosaikCoverage set...we catch the inverse when printing data
    if($cov{$locus}){
	$cov{$locus}{'mapped_reads'} = $1;
	if($cov{$locus}{'mapped_reads'} > 0){
	    $total_mapped_reads = $total_mapped_reads + $cov{$locus}{'mapped_reads'};
	}
    }
    else{
	die "Found a locus $locus in Mosaik asembler file that was not mentioned in the Mosaik coverage file!\n";
    }	

    $line = <AFILE>;
    chomp $line;
}

close(AFILE);

print "\nTOTAL MAPPED READS (assuming only 1 alignment per read): $total_mapped_reads\n";

#calculate the RPKM for each gene
foreach my $locus (keys %cov){
    if($cov{$locus}{'mapped_reads'} > 0){
	$cov{$locus}{'RPKM'} = ($cov{$locus}{'mapped_reads'}/($cov{$locus}{'gene_length'}/1000))/($total_mapped_reads/1000000);
	}
	else {
	    $cov{$locus}{'RPKM'} = 0;
	    $cov{$locus}{'mapped_reads'} = 0;
	}
}

# print values for all genes
&print_output(\%cov);

sub print_output {
	my ($cov) = @_;
	
	open OUT, ">$options{output_file}.txt";
	print OUT "locus\tgene_length\tcoverage_length\tpercentage\tmean_coverage\tRPKM\n";
	
	foreach my $locus (sort {$a cmp $b} keys %$cov) {	    
		print OUT "$locus\t$cov->{$locus}->{'gene_length'}\t$cov->{$locus}->{'coverage_length'}\t$cov->{$locus}->{'percentage'}\t$cov->{$locus}->{'mean_coverage'}\t$cov->{$locus}->{'RPKM'}\n";
	}
	
	close OUT;
}

sub check_parameters {
    my $options = shift;
    
    # make sure required arguments were passed
    my @required = qw( coverage_file output_file assembler_file);
    for my $option ( @required ) {
        unless  ( defined $$options{$option} ) {
            die "--$option is a required option";
        }
    }
}
