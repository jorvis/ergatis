#!/usr/bin/env perl

=head1 NAME

combine_lgt_phylo_para_evidence.pl - Description

=head1 SYNOPSIS

 USAGE: combine_lgt_phylo_para_evidence.pl
       --input_file=/path/to/some/input.file
       --output=/path/to/transterm.file
     [ --log=/path/to/file.log
       --debug=3
       --help
     ]

=head1 OPTIONS

B<--class_c_genes_file,-c>

B<--alien_hunter_file,-a>

B<--sighunt_file,-s>

B<--selfsim_file,-S>

B<--selfsim_cutoff>
  The cutoff chi-square value for a given sliding window

B<--sliding_window_size,-w>
	The size of the sliding window used in selfsim

B<--output_file,-o>

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

my $SELFSIM_CUTOFF = 50;
my $WINDOW_SIZE = 5000; # Applies to selfsim
####################################################

my %options;

# Allow program to run as module for unit testing if necessary (change extension to .pm)
main() unless caller();

sub main {
    my $results = GetOptions (\%options,
						 "class_c_genes_file|c=s",
						 "alien_hunter_file|a=s",
						 "sighunt_file|s=s",
						 "selfsim_file|S=s",
						 "window_size|w=i",
						 "selfsim_cutoff=i",
                         "output_file|o=s",
                         "log|l=s",
                         "debug|d=s",
                         "help|h"
                          );

    &check_options(\%options);

	my $selfsim_cutoff = defined $options{'selfsim_cutoff'} ? $options{'selfsim_cutoff'} : $SELFSIM_CUTOFF;
	my $window_size = defined $options{'window_size'} ? $options{'window_size'} : $WINDOW_SIZE;

	# Parse genes that have been classified as potential class-C LGT
	my $genes = parse_class_c_genes($options{'class_c_genes_file'});

	# If the parametric output file was passed in, parse the file
	parse_alien_hunter_data($options{'alien_hunter_file'}, $genes) if defined $options{'alien_hunter_file'};
	parse_selfsim_data($options{'selfsim_file'}, $selfsim_cutoff, $sliding_window_size, $genes) if defined $options{'selfsim_file'};
	parse_sighunt_data($options{'sighunt_file'}, $genes) if defined $options{'sighunt_file'};

	my $outfh
	open $outfh, ">".$options{'output_file'} or &_log($ERROR, "Cannot open " .$options{'output_file'}. " for writing: $!")
	close $outfh;
}

sub parse_class_c_genes {
	my $file = shift;

	my %genes_hash;
	open IN, $file or &_log($ERROR, "Cannot open $file for reading: $!");
	while (<IN>) {
		chomp;
		my @fields = split;
		# Save start and end coordinates per LGT gene hit
		$genes_hash{$fields[0]}{'start'} = $fields[1];
		$genes_hash{$fields[0]}{'end'} = $fields[2];
	}
	close IN;

	return \%genes_hash;
}

sub parse_alien_hunter_file {
	my $file = shift;
	my $genes_h = shift;

}

sub parse_selfsim_data {
	my $file = shift;
	my $cutoff = shift;
	my $window = shift;
	my $genes_h = shift;

  my $pos_to_shift = $window / 2;

	open IN, $file or &_log($ERROR, "Cannot open $file for reading: $!");
	while (<IN>) {
		chomp;
		next if ($_ =~ /^#/);
		my @fields = split;
		# Only study those where the chi-sq value is less than the cutoff
		if (@fields[1] >= $cutoff) {
			my $window_start = $fields[0] - $pos_to_shift + 1; # want to start at position 1
			my $window_end = $fields[0] + $pos_to_shift;
			foreach my $gene (keys %$genes_h) {
				# Make flag if window completely envelopes gene boundaries
				$genes_h->{$gene}->{'selfsim'} = 1 if ( $window_start <= $genes_h->{$gene}->{'start'} && $genes_h->{$gene}->{'end'} <= $window_end );
			}
		}
	}
	close IN;
}

sub parse_sighunt_data {
	my $file = shift;
	my $genes_h = shift;

	open IN, $file or &_log($ERROR, "Cannot open $file for reading: $!");
	while (<IN>) {
		chomp;
		my @fields = split;
	}
	close IN;
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

    foreach my $req ( qw( class_c_genes_file output_file) ) {
        &_log($ERROR, "Option $req is required") unless( $opts->{$req} );
    }

	&_log($ERROR, "Must pass at least one output from among AlienHunter, Selfsim, or Sighunt") unless (defined $opts->{'alien_hunter_file'} || defined $opts->{'sighunt_file'} || defined $opts->{'selfsim_file'});
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
