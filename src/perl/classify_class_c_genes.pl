#!/usr/bin/env perl

=head1 NAME

classify_class_c_genes.pl - Calculate h-score of two blast files and filter out the class-c genes

=head1 SYNOPSIS

 USAGE: classify_class_c_genes.pl
       --sister_file=/path/to/sister/blast.m8
       --outside_file=/path/to/outside/blast.m8
       --output_dir=/path/to/class_c.tsv
     [ --bitscore_thresh=100
	   --hscore_thresh=30
	   --log=/path/to/file.log
       --debug=3
       --help
     ]

=head1 OPTIONS

B<--sister_file>
	The m8-formatted BLAST output where the subject database was restricted to only the GI accessions in question
	File should only have the best hit for each query gene

B<--outside_file>
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
                         "sister_file=s",
                         "outside_file=s",
                         "output_dir|o=s",
						 "bitscore_thresh:100",
						 "hscore_thresh:30",
                         "log|l=s",
                         "debug|d=s",
                         "help|h"
                          );

    &check_options(\%options);
	parse_m8(\%genes, $options{'sister_file'}, $options{'outside_file'});
	calculate_h_score(\%genes);
	print_class_c_genes(\%genes, \%options);
}

sub parse_m8 {
	my ($gene_h, $only_f, $exclude_f) = @_;
	open SISTER, $only_f || die("Cannot open $only_f for reading: $!");
	while (<SISTER>) {
		chomp;
		my @hit = split(/\t/);
		$gene_h->{$hit[0]}->{'sister_subj'} = $hit[1];
		$gene_h->{$hit[0]}->{'sister_bit'} = $hit[11];
	}
	close SISTER;
	open OUTSIDE, $exclude_f || die("Cannot open $exclude_f for reading: $!");
	while (<OUTSIDE>) {
		chomp;
		my @hit = split(/\t/);
		$gene_h->{$hit[0]}->{'outside_subj'} = $hit[1];
		$gene_h->{$hit[0]}->{'outside_bit'} = $hit[11];
	}
	close OUTSIDE;
}

sub calculate_h_score {
	my $gene_h = shift;
	foreach my $hit (keys %$gene_h) {
		my $only = $gene_h->{$hit}->{'sister_bit'};
		my $exclude = $gene_h->{$hit}->{'outside_bit'};
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
	#my $all_file = $outdir."/all_genes.tsv";
	my $b_thresh = $opts->{'bitscore_thresh'};
	my $h_thresh = $opts->{'hscore_thresh'};
	open CFH, ">".$class_c_file || die("Cannot open $class_c_file for writing: $!");
	#open AFH, ">".$all_file || die("Cannot open $all_file for writing: $!");
	print CFH "query\tsister_hit\toutside_hit\tsister_bit\toutside_bit\th_score\n";
	#print AFH "query\tsister_hit\toutside_hit\tsister_bit\toutside_bit\th_score\thighest_lgt_class\n";
	foreach my $hit (keys %$gene_h) {
		my $gene_str = "$hit\t" . $gene_h->{$hit}->{'sister_subj'} . "\t" . $gene_h->{$hit}->{'outside_subj'} . "\t" . $gene_h->{$hit}->{'sister_bit'} . "\t" . $gene_h->{$hit}->{'outside_bit'} . "\t" . $gene_h->{$hit}->{'h_score'};
		#print AFH $gene_str;

		# If gene meets class C thresholds print to that file, and mark it is a class C gene in "all genes" file
		if ($gene_h->{$hit}->{'h_score'} >= $h_thresh && $gene_h->{$hit}->{'outside_bit'} >= $b_thresh){
			print CFH $gene_str . "\n";
			#print AFH "\tC\n";
		} else {
			#print AFH "\t-\n";
		}
	}
	close CFH;
	#close AFH;
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

    foreach my $req ( qw(sister_file outside_file output_dir) ) {
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
