#!/usr/bin/perl

=head1 NAME

deseq.pl - run the deseq package to compare RNA-seq data

=head1 SYNOPSIS

  USAGE: deseq.pl [
            --r_path=path to R package executable
            --deseq_path=path to the DEseq R script
            --sample_counts=a file containing the replicate id, phenotype, and counts file path for each sample to compare
            --output_dir=output directory
            --annotation=optional file with functional annotation for each gene
            --more_options=additional options
            --help
          ]
=cut

use strict;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Pod::Usage;
use POSIX;

my %options = ();
my $results = GetOptions (\%options, 
			  'r_path=s',
			  'deseq_path=s',
			  'sample_counts=s',
			  'output_dir=s',
			  'annotation=s',
			  'more_options=s',
                          'help|h') || pod2usage();

## display documentation
if( $options{'help'} ){
    pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDERR} );
}

## make sure everything passed was peachy
&check_parameters(\%options);

my $r_bin = $options{'r_path'};
my $deseq_script = $options{'deseq_path'};
my $counts = $options{'sample_counts'};
my $odir = $options{'output_dir'};
my $annotation = $options{'annotation'};
my $more_options = $options{'more_options'};

print "\nRunning $r_bin $more_options --args $counts $odir $annotation < $deseq_script";

system("$r_bin $more_options --args $counts $odir $annotation < $deseq_script");

sub check_parameters {
    my $options = shift;
    
    ## make sure R package path, Deseq script path, counts, output dir were passed
    unless ( $options{r_path} && $options{deseq_path} && $options{sample_counts} && $options{output_dir} ) {
        die "You must pass --r_path --deseq_path --sample_counts and --output_dir arguments!";
    }
    
    ## handle some defaults
    $options{annotation} = "" unless ($options{annotation});
    $options{more_options} = "--slave --vanilla" unless ($options{more_options});
}

