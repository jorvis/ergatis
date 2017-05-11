#!/usr/bin/env perl

=head1 NAME

classify_class_a_genes.pl - Determine if class C genes are class B LGT genes by the avg h-score of their ortholog group

=head1 SYNOPSIS

 USAGE: classify_class_a_genes.pl
       --input_file=/path/to/bitscore/file
       --output_dir=/path/to/dir/
     [ --bitscore_thresh=100
	   --log=/path/to/file.log
       --debug=3
       --help
     ]

=head1 OPTIONS

B<--input_file,-i>
	File that lists bitscores of all genes and B-class labelings from classify_class_b_genes.pl

B<--output_dir,-o>
	Directory to write output into

B<--bitscore_thresh>
	Keep genes whose "only" bitscore is under this value.  Default is 100

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
                         "output_dir|o=s",
						 "bitscore_thresh:100",
                         "log|l=s",
                         "debug|d=s",
                         "help|h"
                          );

    &check_options(\%options);
	parse_input_file(\%genes, $options{'input_file'});
	print_class_a_genes(\%genes, \%options);
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
		$gene_h->{$hit[0]}->{'cluster'} = $hit[6];
		$gene_h->{$hit[0]}->{'ortholog_h_score'} = $hit[7];
		$gene_h->{$hit[0]}->{'lgt_class'} = $hit[8];
	}
	close IFH;
}

sub print_class_a_genes {
	my ($gene_h, $opts) = @_;
	my $outdir = $opts->{'output_dir'};
	my $class_a_file = $outdir."/class_a.tsv";
	my $all_file = $outdir."/all_genes.tsv";
	my $bit_thresh = $opts->{'bitscore_thresh'};
	open AFH, ">".$class_a_file || &_log($ERROR, "Cannot open $class_a_file for writing: $!");
	open FH, ">".$all_file || &_log($ERROR, "Cannot open $all_file for writing: $!");
	print AFH "query\tonly_hit\texclude_hit\tonly_bit\texclude_bit\th_score\tortholog_cluster\tortholog_h_score\n";
	print FH "query\tonly_hit\texclude_hit\tonly_bit\texclude_bit\th_score\tortholog_cluster\tortholog_h_score\thighest_lgt_class\n";
	foreach my $hit (keys %$gene_h) {
		my $gene_str = "$hit\t" . $gene_h->{$hit}->{'only_subj'} . "\t" . $gene_h->{$hit}->{'exclude_subj'} . "\t" . $gene_h->{$hit}->{'only_bit'} . "\t" . $gene_h->{$hit}->{'exclude_bit'} . "\t" . $gene_h->{$hit}->{'h_score'} . "\t" . $gene_h->{$hit}->{'cluster'} . "\t" . $gene_h->{$hit}->{'ortho_h'};
		print AFH $gene_str;

		my $lgt_class = $gene_h->{$hit}->{'lgt_class'};
		# If class B gene has "only_bit" < 100 and all ortholog members have "only_bit" < 100, then it is a class A gene
		if ($lgt_class eq 'B' && $gene_h->{$hit}->{'only_bit'} < $bit_thresh){
			my $class_a_flag = 1;
			# Verify all other members in this cluster have 'only' bitscores under the threshold
			foreach my $ortho_hit (keys %$gene_h) {
				if ($gene_h->{$hit}->{'cluster'} eq $gene_h->{$hit}->{'cluster'} && $gene_h->{$ortho_hit}->{'only_bit'} >= $bit_thresh) {
					$class_a_flag = 0;
					last
				}
			}
			if ($class_a_flag) {
				$lgt_class = "A";
				print AFH $gene_str . "\n";
			}
		}
		print FH "\t$lgt_class\n";
	}
	close AFH;
	close FH;
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
