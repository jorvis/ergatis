#!/usr/bin/env perl
=head1 NAME

filter_dups_lc_seqs.pl - Run Picard tools & Prinseq to filter out duplicate and low-complexity sequences

=head1 SYNOPSIS

 USAGE: filter_dups_lc_seqs.pl
       --input_file=/path/to/some/input.bam
       --output_dir=/path/to/output/dir
	   --prinseq_path=/path/to/prinseq.pl
	   --samtools_path=/path/to/samtools/exec
	   --picard_path=/path/to/picard/dir
     [ --prefix='prefix'
	   --rm_duplicates
	   --rm_low_complexity
       --lc_method='dust'
       --lc_threshold=7
       --tmp_dir/path/to/scratch
	   --log=/path/to/file.log
       --debug=3
       --help
     ]

=head1 OPTIONS

B<--input_file,-i>
	A paired-end BAM file

B<--output_dir,-o>

B<--tmp_dir,-t>

B<--prefix,-p>
	Abbreviation used to name post-processing tab file

B<--prinseq_path,-P>
	The path to the prinseq-lite Perl script

B<--samtools_path,-s>
	The path to the SAMtools executable

B<--picard_path,-P>
	The path to the Picard JAR file

B<--rm_duplicates,-d>
	If flag is present, will remove duplications of sequences

B<--rm_low_complexity,-c>
	If flag is present, will remove sequences deemed to have low complexity

B<--lc_method,-m>
    Method used to filter out low complexity seqs.  Default is 'dust'

B<--lc_threshold,-t>
    Threshold for filtering low complexity seqs.  Default is 7.

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
use LGT::LGTSeek;

############# GLOBALS AND CONSTANTS ################
my $debug = 1;
my ($ERROR, $WARN, $DEBUG) = (1,2,3);
my $logfh;
my $DEFAULT_LC_THRESHOLD = 7;
my $DEFAULT_LC_METHOD = "dust";
my $JAVA_PATH = "/usr/local/packages/java-jre-1.8.0/bin/java";	# Picard-2.0.1 needs Java 1.8
#my $JAVA_OPTS = '-Xmx2g';
####################################################

my %options;
my $results = GetOptions (\%options,
                         "input_file|i=s",
                         "output_dir|o=s",
                         "tmp_dir|t=s",
						 "prefix|p=s",
						 "prinseq_path|P=s",
						 "samtools_path|s=s",
						 "picard_path|p=s",
						 "rm_duplicates|d",
						 "rm_low_complexity|c",
                         "lc_method|m=s",
                         "lc_threshold|t=i",
                         "log|l=s",
                         "debug=s",
                         "help|h"
                          );

&check_options(\%options);

my $count = check_empty_file();
if ($count == 0) {
	print STDERR "Input file $options{'input_file'} has no alignments.  Exiting and not creating output.\n";
	exit(0);
}

my $tmp_dir = (defined $options{'tmp_dir'}) ? $options{'tmp_dir'} : $options{'output_dir'} . "/tmp";

my $lc_method = defined $options{'lc_method'} ? $options{'lc_method'} : $DEFAULT_LC_METHOD;
my $lc_threshold = defined $options{'lc_threshold'} ? $options{'lc_threshold'} : $DEFAULT_LC_THRESHOLD;
my $prefix = (defined $options{'prefix'}) ? $options{'prefix'} : "post_prinseq";

my $lgtseek = LGT::LGTSeek->new( {
		'prinseq_bin'	=> $options{'prinseq_path'},
		'Picard_jar'	=> $options{'picard_path'},
		'java_bin'		=> $JAVA_PATH,
		'java_opts'		=> '',
		'samtools_bin'	=> $options{'samtools_path'},
		'paired_end'	=> 1,
		'verbose'		=> 1
	} );

my $filtered_bam = $lgtseek->prinseqFilterBam( {
		'input_bam' 	=>	$options{'input_file'},
		'tmp_dir'		=>	$tmp_dir,
		'output_dir'	=>	$options{'output_dir'},
		'lc_method'		=>	$lc_method,
		'lc_threshold'	=>	$lc_threshold,
		'dedup'			=>	$options{'rm_duplicates'},
		'rm_low_cmplx'	=>	$options{'rm_low_complexity'}
	} );

my $counts_file = $options{'output_dir'} . "/" . $prefix . ".post_processing.tab";
my (@header, @vals);
push( @header, 'lgt_pass_prinseq' );
push( @vals,   $filtered_bam->{count} );

LGT::LGTSeek->print_tab( $counts_file, \@header, \@vals );
exit(0);

# Check if input BAM file is empty and return if any output is in head
# The chk_empty subroutine in LGT::LGTSeek works but we want to silently exit without making output
sub check_empty_file {
	return `$options{samtools_bin} view $options{'input_file'} | head | wc -l`;
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

   foreach my $req ( qw(input_file output_dir prinseq_path samtools_path picard_path) ) {
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
