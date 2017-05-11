#!/usr/bin/env perl

=head1 NAME

classify_class_b_genes.pl - Determine if class C genes are class B LGT genes by the avg h-score of their ortholog group

=head1 SYNOPSIS

 USAGE: classify_class_b_genes.pl
       --input_file=/path/to/hscore/file
	   --orthlogs_file=/path/to/ortholog/clusters.txt
       --output_dir=/path/to/dir/
     [ --hscore_thresh=30
	   --log=/path/to/file.log
       --debug=3
       --help
     ]

=head1 OPTIONS

B<--input_file,-i>
	File that lists h_scores of all genes and C-class labelings from classify_class_c_genes.pl

B<--orthologs_file,-O>
	File from FastOrtho that lists the end-result ortholog clusters

B<--output_dir,-o>
	Directory to write output into

B<--hscore_thresh>
	Keep genes whose ortholog group's average h-score exceeds this value.  Default is 30

B<--log,-l>
    Logfile.

B<--debug,-d>
    1,2 or 3. Higher values more verbose.

B<--help,-h>
    Print this message

=head1  DESCRIPTION

 DESCRIPTION
 
=head1  INPUT

    Describe the input

=head1 OUTPUT

    Describe the output

=head1  CONTACT

    Shaun Adkins
    sadkins@som.umaryland.edu

=cut

use strict;
use warnings;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Pod::Usage;

############# GLOBALS AND CONSTANTS ################
my $debug = 1;
my ($ERROR, $WARN, $DEBUG) = (1,2,3);
my $logfh;
####################################################

my %options;
my %genes;
my %ortho;

# Allow program to run as module for unit testing if necessary (change extension to .pm)
main() unless caller();
exit(0);

sub main {
    my $results = GetOptions (\%options,
						 "input_file|i=s",
						 "orthologs_list|O=s",
                         "output_dir|o=s",
						 "hscore_thresh:30",
                         "log|l=s",
                         "debug|d=s",
                         "help|h"
                          );

    &check_options(\%options);
	parse_input_file(\%genes, $options{'input_file'});
	parse_orthologs(\%ortho, $options{'orthologs_file'});
	calculate_avg_ortho_h_score(\%genes, \%ortho);
	print_class_b_genes(\%genes, \%options);
}

sub parse_input_file {
	my ($gene_h, $input) = @_;
	open IFH, $input || &_log($ERROR, "Cannot open $input for reading: $!");
	while (<IFH>) {
		chomp;
		my @hit = split(/\t/);
		$gene_h->{$hit[0]}->{'only_subj'} = $hit[1];
		$gene_h->{$hit[0]}->{'exclude_subj'} = $hit[2];
		$gene_h->{$hit[0]}->{'only_bit'} = $hit[3];
		$gene_h->{$hit[0]}->{'exclude_bit'} = $hit[4];
		$gene_h->{$hit[0]}->{'h_score'} = $hit[5];
		$gene_h->{$hit[0]}->{'lgt_class'} = $hit[6];
	}
	close IFH;
}

sub parse_orthologs {
	my ($ortho_h, $ortho_file) = @_;
	open ORTHO, $ortho_file || &_log($ERROR, "Cannot open $ortho_file for reading: $!");
	while (my $line = <ORTHO>) {
		chomp $line;
		# Parse out ortholog ID and genes
		if ($line =~ /^(\S+)\s\(\d+\sgenes\d+\staxa\):\s+(.+)$/) {
			my $cluster_id  = $1;
			my $gene_str = $2;
			my @gene_arr = split(/\s+/, $gene_str);
			$ortho_h->{$cluster_id} = \@gene_arr;
		} else {
			&_log($WARN, "Line $line did not match the regex.  Please look into this");
		}
	}
	close ORTHO;
}

sub calculate_avg_ortho_h_score {
	my ($genes_h, $ortho_h) = @_;

	foreach my $ortholog (keys %$ortho_h) {
		my $h_score_sum = 0;
		my $num_genes = scalar @{$ortho_h->{$ortholog}};
		foreach my $gene (@{$ortho_h->{$ortholog}}){
			# Get rid of species name
			$gene =~ s/\(.*$//g;
			$h_score_sum += $genes_h->{$gene}->{'h_score'};
		}
		my $h_score = $h_score_sum/$num_genes;
		# Insert key for average ortholog h_score for each gene in that ortholog group (for easy retrieval later)
		foreach my $gene (@{$ortho_h->{$ortholog}}){
			$genes_h->{$gene}->{'cluster'} = $ortholog;
			$genes_h->{$gene}->{'ortho_h'} = $h_score;
		}
	}
}

sub print_class_b_genes {
	my ($gene_h, $opts) = @_;
	my $outdir = $opts->{'output_dir'};
	my $class_b_file = $outdir."/class_b.tsv";
	my $all_file = $outdir."/all_genes_w_hscore.tsv";
	my $h_thresh = $opts->{'hscore_thresh'};
	open BFH, ">".$class_b_file || &_log($ERROR, "Cannot open $class_b_file for writing: $!");
	open AFH, ">".$all_file || &_log($ERROR, "Cannot open $all_file for writing: $!");
	print BFH "query\tonly_hit\texclude_hit\tonly_bit\texclude_bit\th_score\tortholog_cluster\tortholog_h_score\n";
	print AFH "query\tonly_hit\texclude_hit\tonly_bit\texclude_bit\th_score\tortholog_cluster\tortholog_h_score\thighest_lgt_class\n";
	foreach my $hit (keys %$gene_h) {
		my $gene_str = "$hit\t" . $gene_h->{$hit}->{'only_subj'} . "\t" . $gene_h->{$hit}->{'exclude_subj'} . "\t" . $gene_h->{$hit}->{'only_bit'} . "\t" . $gene_h->{$hit}->{'exclude_bit'} . "\t" . $gene_h->{$hit}->{'h_score'} . "\t" . $gene_h->{$hit}->{'cluster'} . "\t" . $gene_h->{$hit}->{'ortho_h'};
		print AFH $gene_str;

		my $lgt_class = $gene_h->{$hit}->{'lgt_class'};
		# If gene meets class C thresholds print to that file, and mark it is a class C gene in "all genes" file
		if ($lgt_class eq 'C' && $gene_h->{$hit}->{'ortho_h'} >= $h_thresh){
			$lgt_class = "B";
			print BFH $gene_str . "\n";
		}
		print AFH "\t$lgt_class\n";
	}
	close BFH;
	close AFH;
}

sub check_options {
    my $opts = shift;
    if( $opts->{'help'} ) {
        &_pod;
    }

    if( $opts->{'log'} ) {
        open( $logfh, "> $opts->{'log'}") or die("Can't open log file ($!)");
    }

    $debug = $opts->{'debug'} if( $opts->{'debug'} );

    foreach my $req ( qw(only_file exclude_file output_dir) ) {
        &_log($ERROR, "Option $req is required") unless( $opts->{$req} );
    }
}

sub _log {
    my ($level, $msg) = @_;
    if( $level <= $debug ) {
        print STDOUT "$msg\n";
    }
    print $logfh "$msg\n" if( defined( $logfh ) );
    exit(1) if( $level == $ERROR );
}

sub _pod {
    pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDERR} );
}
