#!/usr/bin/env perl
=head1 NAME

filter_low_scoring_bers.pl - Filters out BER entries with low-scoring identity and similarity

=head1 SYNOPSIS

 USAGE: filter_low_scoring_bers.pl
       --input_file=/path/to/some/ber.nr
	   --btab_file=/path/to/some/ber.nr.btab
       --output=/path/to/filtered_ber.nr
	   --identity=0.35
	   --similarity=0.55
     [ --log=/path/to/file.log
       --debug=3
       --help
     ]

=head1 OPTIONS

B<--input_file,-i>
	BER results file

B<--btab_file, -b>
	Corresponding BER btab file

B<--output_file,-o>
	Path to new BER results file.  A .btab extension will be applied to create the new .btab file also

B<--identity, -i>
	Cutoff to filter percent identity.  Range from 0.0 to 1.0

B<--similiarity, -s>
	Cutoff to filter percent similarity. Range from 0.0 to 1.0

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

my $identity = 0.35;
my $similarity = 0.55;
####################################################

my ($ber_file, $btab_file, $out_file, $out_btab_file);
my %options;

my $results = GetOptions (\%options,
                "input_file|f=s",
				"btab_file|b=s",
                "output_file|o=s",
				"identity|i=i",
				"similarity|s=i",
                "log|l=s",
                "debug|d=s",
                "help|h"
                );

&check_options(\%options);
process_ber_file($ber_file);

# Parse the whole BER output file
sub process_ber_file {
	my ($in_file) = @_;
	my @query_arr;

	open IN, $in_file or die "Cannot open $in_file for reading: $!\n";
	open my $out_btab_fh, ">$out_btab_file" or die "Cannot open $out_btab_file for writing: $!\n";
	open my $outfh, ">$out_file" or die "Cannot open output file $out_file for writing: $!\n";
	while (<IN>) {
		my $line = $_;
		chomp $line;
		if ($line =~ /^Query = (.+)/) {
			if (scalar @query_arr > 0) {
				process_query(\@query_arr, $out_btab_fh, $outfh);
				@query_arr = ();
			}
		}
		push @query_arr, $line;
	}
	process_query(\@query_arr, $out_btab_fh, $outfh);

	close $outfh;
	close $out_btab_fh;
	close IN;
}

# Process all the data related to a particular query ID
sub process_query {
	my $q = shift;
	my $out_btab_fh = shift;
	my $outfh = shift;

	my @kept_hits;	# array for kept subject IDs
	my $query;
	my @q_arr = @$q;

	for my $i (0..$#q_arr) {
		my $line = $q_arr[$i];
		if ($line =~ /^Query = (.+)/) {
			$query = $1;
			print $outfh "$line\n";
		} elsif ($line =~ /^(\w+)\s+(\d+\.\d+)\s+(\d+\.\d+)\s+(\d+\.\d+)$/) {
            # Matching Subject_id   %Match  %Id %Sim
            my $subj = $1; 
            my $id = $3; 
            my $sim = $4; 
            if ($id >= $identity && $sim >= $similarity) {
				push @kept_hits, $subj;
				print $outfh "$line\n";
                my $btab_line = grab_btab_line($btab_file, $query, $subj);
                print $out_btab_fh "$btab_line\n";
            }
        } elsif ($line =~ /^$query\(/) {
			_log($DEBUG, "Found " . scalar(@kept_hits) . " number of good hits for query $query");
			# When we encounter first alignment section of query, break off and process separately
			my @alignment_arr = splice(@q_arr, $i);
			process_alignments(\@alignment_arr, \@kept_hits, $query, $out_btab_fh, $outfh);
			return;
		} else {
			# All other misc. lines should be printed
			print $outfh "$line\n";
		}
	}
}

# Process the alignment portion of a given query sequence
sub process_alignments {
	my ($a, $kept, $query, $out_btab_fh, $outfh) = @_;
	my @single_alignment;	
	my $good_subj = 0;
	my @a_arr = @$a;

	foreach my $i (0..$#a_arr){
		my $line = $a_arr[$i];
		# Break into signle query/subject alignments, and check if subject was filtered
		if ($line =~ /^$query\(/) {
			my $next_line = $a_arr[$i+1];
			foreach my $hit (@$kept) {
				if ($next_line =~ /^$hit\(/) {
					$good_subj = 1;
					_log($DEBUG, "Found good alignment for $query: $hit");
					last;
				}
			}
			if (scalar @single_alignment > 0) {
				# Only print prev alignment if this was a good filtered hit
				if ($good_subj) {
					print $outfh join("\n", @single_alignment) . "\n";	
				}
				@single_alignment = ();
				$good_subj = 0;
			}
		}
		push @single_alignment, $line;
	}

    if ($good_subj) {
        print $outfh join("\n", @single_alignment) . "\n";
	}

}

sub grab_btab_line {
	my ($file, $q, $s) = @_;
	open BTAB, $file or die "Cannot open btab file $file for reading: $!\n";
	while (<BTAB>){
		my $line = $_;
		chomp $line;
		# A bit iffy using ".+" in the regex but it allows for usage of the lookbehind
		return $line if ($line =~ /^$q\s+.+(?<=minidb)\s+$s\s+/);
	}
	close BTAB;
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

   foreach my $req ( qw(input_file btab_file output_file) ) {
       &_log($ERROR, "Option $req is required") unless( $opts->{$req} );
   }

   $ber_file = $opts->{'input_file'};
   $btab_file = $opts->{'btab_file'};
   $out_file = $opts->{'output_file'};
   $out_btab_file = $out_file . ".btab";
   $identity = $opts->{'identity'} if $opts->{'identity'};
   $similarity = $opts->{'similarity'} if $opts->{'similarity'};
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
