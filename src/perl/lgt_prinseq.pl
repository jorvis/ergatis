#!/usr/bin/env perl
package run_prinseq;    # Change this to reflect name of script.  If you want to use as module, make sure to save as .pm instead of .pl

=head1 NAME

run_prinseq.pl - Run Prinseq to filter out duplicate and low-complexity sequences

=head1 SYNOPSIS

 USAGE: perl_template.pl
       --input_file=/path/to/some/input.file
       --output=/path/to/transterm.file
	   --prinseq_path=/path/to/prinseq.pl
	   --samtools_path=/path/to/samtools/exec
	   --picard_path=/path/to/picard/dir
     [ --rm_duplicates
	   --rm_low_complexity
	   --log=/path/to/file.log
       --debug=3
       --help
     ]

=head1 OPTIONS

B<--input_file,-i>

B<--output_file,-o>

B<--prinseq_path,-p>
	The path to the prinseq-lite Perl script

B<--samtools_path,-s>
	The path to the SAMtools executable

B<--picard_path,-P>
	The path to the Picard java suite

B<--rm_duplicates,-d>
	If flag is present, will remove duplications of sequences

B<--rm_low_complexity,-c>
	If flag is present, will remove sequences deemed to have low complexity

B<--log,-l>
    Logfile.

B<--debug>
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

# Allow program to run as module for unit testing if necessary
main() unless caller();

sub main {
	my $results = GetOptions (\%options,
                         "input_file|i=s",
                         "output_file|o=s",
						 "prinseq_path|p=s",
						 "samtools_path|s=s",
						 "picard_path|p=s",
						 "rm_duplicates|d",
						 "rm_low_complexity|c",
                         "log|l=s",
                         "debug=s",
                         "help|h"
                          );

    &check_options(\%options);

	# Filter BAM input by IDs
	# Generate concatenated fastq files in order to remove duplicates
	# Generate single-read fastq files for low_complexity filtering
	# Merge IDs obtained from filtering dups and low_complexity
	# Lastly, filter out identified duplicates and low-complexity seqs
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

   foreach my $req ( qw(input_file output_file prinseq_path samtools_path picard_path) ) {
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
