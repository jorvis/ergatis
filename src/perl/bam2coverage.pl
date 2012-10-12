#!/usr/bin/perl

eval 'exec /usr/bin/perl  -S $0 ${1+"$@"}'
    if 0; # not running under some shell

=head1 NAME

bam2coverage.pl - Generate coverage stats for alignments in BAM files.  
                - If given a gff3 file, also computes expression of genes. 
	            - Produces a chromosomes.expr,genes.expr, and genes.counts file.

=head1 SYNOPSIS

 USAGE: bam2coverage.pl
	--bam_file=/path/to/bam/file
	--total_mapped_reads= optionally pass in precalculated total mapped reads for rpkm calculation
	--stranded= 1 or 0 to indicate the library was stranded (default: 0)
        --reverse_strand = 1 or 0, to reverse the strand of features, useful for wierd RNA-seq protocols [default 0]
	--gff3_file=/path/to/gff3/file
	--feature_type=feature to use for transcripts (default: gene)
	--genome_only= 1 or 0 calculate only chromosome level coverage information
	--features_only= 1 or 0  calculate only feature level coverage information
	--output_stub=path and filename stub to use for output files
	[ --log=/path/to/file.log
	--help
	]

=head1  CONTACT

	Todd Creasy
    tcreasy@som.umaryland.edu

=cut

use strict;
use warnings;
use Data::Dumper;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use FindBin;
use lib "$FindBin::Bin";
use Pod::Usage;
use NGS::GFF3;
use NGS::Alignment;

my %options = ();
my $results = GetOptions (\%options, 
			  'bam_file=s',
			  'total_mapped_reads=s',
			  'stranded=s',
			  'reverse_strand=s',
			  'gff3_file=s',
			  'feature_type=s',
                          'genome_only=s',
			  'features_only=s',
			  'output_stub=s',
                          'log=s',
                          'help|h') || pod2usage();

# display documentation
if( $options{'help'} ) {
    pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDERR} );
}

# make sure parameters are correct
check_parameters(\%options);

# default is to use NON strand-specific data
if( !defined $options{stranded} ) {
	$options{stranded} = 0;
}

if( !defined $options{reverse_strand} ) {
	$options{reverse_strand} = 0;
}

# get coverage for genome only
if( !defined $options{genome_only} ) {
	$options{genome_only} = 0;
}

# get coverage for specified features only
if( !defined $options{features_only} ) {
	$options{features_only} = 0;
}

# default feature type is gene
if( !defined $options{feature_type} ) {
	$options{feature_type} = "gene";
}

# set default output file name prefix to use --bam_file stub
if( !defined $options{output_stub} ) {
	# regex to pull out stub from file name
	# Example: /my/directory/37_default/test.sorted.bam = test
	($options{output_stub}) = ($options{bam_file} =~ /.*\/([^\.]*)\.*/);
}

# open the log if requested
my $logfh;
if (defined $options{log}) {
    open($logfh, ">".$options{output_stub}.$options{log}) || die "can't create log file: $!";
}

my $calculate_genome_cov = 0;
my $calculate_feature_cov = 0;

if( $options{genome_only} || (!$options{genome_only} && !$options{features_only}) ) {
    $calculate_genome_cov = 1;
}

if( $options{features_only} || (!$options{features_only} && !$options{genome_only}) ) {
    $calculate_feature_cov = 1;
}

if( $calculate_feature_cov && !(defined $options{gff3_file}) ) {
    die "Please provide a GFF3 file with --gff3_file to calculate feature level coverage!";
}

print "\nCalculating coverage from ".$options{bam_file}." stranded = ".$options{stranded}." feature type = ".$options{feature_type}."\n";

# create Alignment object
print "Create alignment object...\n";
my $bam = (defined $options{total_mapped_reads}) ? 
    NGS::Alignment->new( file_type => "BAM", file => $options{bam_file}, total_mapped_reads => $options{total_mapped_reads}) :
    NGS::Alignment->new( file_type => "BAM", file => $options{bam_file} );
print "Done.\n\n";


print "Get total mapped reads...\n";
$options{total_mapped_reads} = $bam->get_total_mapped_reads();
print "Done.\n\n";


print "\nTotal mapped reads: " . $options{total_mapped_reads} . "\n";


if( $calculate_genome_cov ) {

    print "\nCalculating genome coverage...\n";
    my $total_genome_cov = 0;
    my $total_genome_bps = 0;
    my $total_genome_bps_covered = 0;
	
    my $genome_coverage = $bam->get_genome_coverage();
	print "Done.\n\n";

    my @clabels = ("chromosome", "length", "% coverage", "depth", "read_count");
	
    open CCOV, ">".$options{output_stub}.".chromosomes.expr";
    print CCOV join("\t", @clabels, "\n");
	
    foreach my $chrom (@$genome_coverage) {  
		my @data = ($chrom->{'name'},$chrom->{'bps'}, $chrom->{'coverage'}."%", $chrom->{'depth'}."x", $chrom->{'hits'});
		
		$total_genome_cov += $chrom->{'depth'};
		$total_genome_bps += $chrom->{'bps'};
		$total_genome_bps_covered += $chrom->{'bps_covered'};
		
		print CCOV join("\t", @data, "\n");
    }
    close CCOV;

    # output
    print "Total mapped reads: ".$options{total_mapped_reads}."\n";
    print "Total genome length: ".$total_genome_bps."\n";
    print "Total bps covered: ".$total_genome_bps_covered."\n";
    print "Percentage of genome covered: ".(($total_genome_bps_covered/$total_genome_bps)*100)."%\n";
    print "Average coverage depth: ".($total_genome_cov/@$genome_coverage)."x\n";

}


my $features_processed = 0;

if( $calculate_feature_cov ) {
    
    print "\nCalculating feature coverage...\n";

	my $total_mapped_reads = 0;
    my $total_gene_cov = 0;
    my $total_gene_bps = 0;
    my $total_gene_bps_covered = 0;    

    my $genes_coverage;

    # create GFF3 object, to use in our transcripts calculations
    my $gff3 = NGS::GFF3->new( file => $options{gff3_file} );


    #reverse strand on features
    if($options{reverse_strand}){
		print "\nReversing strand on features...";
		$gff3->reverse_feature_strands();
    }
	
	
    if( $options{stranded} ) {
		$genes_coverage = $bam->get_feature_coverages({feature => $options{feature_type}, gff3 => $gff3, stranded => 1});
    } else {
		$genes_coverage = $bam->get_feature_coverages({feature => $options{feature_type}, gff3 => $gff3});
    }
	
    open GCOUNT, ">" . $options{output_stub} . ".genes.counts";
    open GCOV, ">" .  $options{output_stub} . ".genes.expr";
    print GCOV join("\t", "feature id", "length", "hits", "% coverage", "depth", "read_count", "RPKM\n");
	
    foreach my $feature (@$genes_coverage){

		$total_mapped_reads += $feature->{'hits'};
        $total_gene_cov += $feature->{'depth'};
        $total_gene_bps += $feature->{'bps'};
        $total_gene_bps_covered += $feature->{'bps_covered'};
		
        my $rpkm = ($feature->{'bps'} > 0) ? ($feature->{'hits'}/($feature->{'bps'}/1000))/($options{total_mapped_reads}/1000000) : 0;
       
        print GCOV join("\t", 
						$feature->{'name'}, 
						$feature->{'bps'}, 
						$feature->{'hits'}, 
						$feature->{'coverage'}."%", 
						$feature->{'depth'}."x", 
						$feature->{'hits'},
						$rpkm
						)."\n";
		
        print GCOUNT join ("\t", $feature->{'name'}, $feature->{'hits'})."\n";
		
        $features_processed++;		
    }
	
    close GCOV;
    close GCOUNT;
	print "Done.\n\n";

	# output
    print "Processed $features_processed features\n\n";
	
    if( $features_processed > 0 ) {
		print "Total mapped reads: " . $total_mapped_reads . "\n";
		print "Total features length length: ".$total_gene_bps."\n";
		print "Total features bps covered: ".$total_gene_bps_covered."\n";
		print "Percentage of features covered: ".(($total_gene_bps_covered/$total_gene_bps)*100)."%\n";
		print "Average coverage depth: ".($total_gene_cov/$features_processed)."x\n";
    }

}


##############################################

sub _log {
    my $msg = shift;
    print "$msg\n";
    print $logfh "$msg\n" if $logfh;
}

sub check_parameters {
    my $options = shift;
    
    # make sure required arguments were passed
    my @required = qw( bam_file );
    for my $option ( @required ) {
        unless  ( defined $$options{$option} ) {
            die "--$option is a required option";
        }
    }    

	# can only do one type coverage analysis (genome or features)
	if( ((defined $options{genome_only}) && ($options{genome_only} == 1)) && ((defined $options{features_only}) && ($options{features_only} == 1)) ) {
		die "Cannot define genome_only and features_only simultaneously!";
	}
		
}
