#!/usr/bin/env perl

=head1 NAME

classify_class_c_genes.pl - Calculate h-score of two blast files and filter out the class-c genes

=head1 SYNOPSIS

 USAGE: classify_class_c_genes.pl
       --only_file=/path/to/only/blast.m8
       --exclude_file=/path/to/restricted/blast.m8
       --output_dir=/path/to/class_c.tsv
     [ --bitscore_thresh=100
	   --hscore_thresh=30
	   --log=/path/to/file.log
       --debug=3
       --help
     ]

=head1 OPTIONS

B<--only_file>
	The m8-formatted BLAST output where the subject database was restricted to only the GI accessions in question
	File should only have the best hit for each query gene

B<--exclude_file>
	The m8-formatted BLAST output where the subject database excluded those GI accessions from --only_file
	File should only have the best hit for each query gene

B<--output_dir,-o>
	Directory to write output into

B<--bitscore_thresh>
	Keep genes that exceed this bitscore in the --exclude_file.  Default is 100

B<--hscore_thresh>
	Keep genes whose h-score exceeds this value.  Default is 30

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

# Allow program to run as module for unit testing if necessary (change extension to .pm)
main() unless caller();
exit(0);

sub main {
    my $results = GetOptions (\%options,
                         "only_file=s",
                         "exclude_file=s",
                         "output_dir|o=s",
						 "bitscore_thresh:100",
						 "hscore_thresh:30",
                         "log|l=s",
                         "debug|d=s",
                         "help|h"
                          );

    &check_options(\%options);
	parse_m8(\%genes, $options{'only_file'}, $options{'exclude_file'});
	calculate_h_score(\%genes);
	print_class_c_genes(\%genes, \%options);
}

sub parse_m8 {
	my ($gene_h, $only_f, $exclude_f) = @_;
	open ONLY, $only_f || die("Cannot open $only_f for reading: $!");
	while (<ONLY>) {
		chomp;
		my @hit = split(/\t/);
		$gene_h->{$hit[0]}->{'only_subj'} = $hit[1];
		$gene_h->{$hit[0]}->{'only_bit'} = $hit[11];
	}
	close ONLY;
	open EXCLUDE, $exclude_f || die("Cannot open $exclude_f for reading: $!");
	while (<EXCLUDE>) {
		chomp;
		my @hit = split(/\t/);
		$gene_h->{$hit[0]}->{'exclude_subj'} = $hit[1];
		$gene_h->{$hit[0]}->{'exclude_bit'} = $hit[11];
	}
	close EXCLUDE;
}

sub calculate_h_score {
	my $gene_h = shift;
	foreach my $hit (keys %$gene_h) {
		my $only = $gene_h->{$hit}->{'only_bit'};
		my $exclude = $gene_h->{$hit}->{'exclude_bit'};
		if (defined $only && defined $exclude) {
			$gene_h->{$hit}->{'h_score'} = $exclude - $only;
		} else {
			# If either bitscore is missing, print warning and delete hit
			&_log($WARN, "WARNING - Gene $hit was missing a top bitscore from one or both BLAST searches.  Cannot calculate 'h_score'... skipping");
			delete $gene_h->{$hit};
		}
	}
}

sub print_class_c_genes {
	my ($gene_h, $opts) = @_;
	my $outdir = $opts->{'output_dir'};
	my $class_c_file = $outdir."/class_c.tsv";
	my $all_file = $outdir."/all_genes.tsv";
	my $b_thresh = $opts->{'bitscore_thresh'};
	my $h_thresh = $opts->{'hscore_thresh'};
	open CFH, ">".$class_c_file || die("Cannot open $class_c_file for writing: $!");
	open AFH, ">".$all_file || die("Cannot open $all_file for writing: $!");
	print CFH "query\tonly_hit\texclude_hit\tonly_bit\texclude_bit\th_score\n";
	print AFH "query\tonly_hit\texclude_hit\tonly_bit\texclude_bit\th_score\thighest_lgt_class\n";
	foreach my $hit (keys %$gene_h) {
		my $gene_str = "$hit\t" . $gene_h->{$hit}->{'only_subj'} . "\t" . $gene_h->{$hit}->{'exclude_subj'} . "\t" . $gene_h->{$hit}->{'only_bit'} . "\t" . $gene_h->{$hit}->{'exclude_bit'} . "\t" . $gene_h->{$hit}->{'h_score'};
		print AFH $gene_str;

		# If gene meets class C thresholds print to that file, and mark it is a class C gene in "all genes" file
		if ($gene_h->{$hit}->{'h_score'} >= $h_thresh && $gene_h->{$hit}->{'exclude_bit'} >= $b_thresh){
			print CFH $gene_str . "\n";
			print AFH "\tC\n";
		} else {
			print AFH "\t-\n";
		}
	}
	close CFH;
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
