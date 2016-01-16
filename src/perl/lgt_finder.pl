#!/usr/bin/env perl
=head1 NAME

lgt_finder.pl - Wrapper script to find LGT hits between donor and host

=head1 SYNOPSIS

 USAGE: lgt_finder.pl
       --input_file/input_list=/path/to/some/input
       --output_file=/path/to/transterm.txt
	   --ref_lineage="Drosophila"
	 [
	   --output_prefix=blah
	   --max_overlap=20
	   --min_length=0
	   --log=/path/to/file.log
       --debug=3
       --help
     ]

=head1 OPTIONS

B<--input_file,-i>
	A file consisting of the best blast hits per query, obtainable from lgt_best_blast_hit.pl

B<--input_list,-I>
	A list file consisting of input files

B<--output_dir,-o>
	Path name to output directory. 

B<--output_prefix,-p>
	Prefix name to give output files

B<--max_overlap,-M>
	Overlap threshold that too good hits will not exceed.  Default is 20 bases

B<--min_length,-m>
	Threshold for minimum alignment length of a hit.  Default is 0 bases.

B<--ref_lineage,-d>
	Name of the lineage of the reference genome.  This will most likely be the host/recipient genome

B<--log,-l>
    Logfile.

B<--debug>
    1,2 or 3. Higher values more verbose.

B<--help>
    Print this message

=head1  DESCRIPTION

 DESCRIPTION
 
=head1  INPUT

    BlastN results that have been formatted into the m8 format

=head1 OUTPUT

    Three output files containing information about the best hit and LCA for a given query.
	One output will correspond to traces mapping to the host genome
	One output will correspond to traces mapping to the donor genome
	The final output will correspond to the overall mate clones

=head1  CONTACT

    Shaun Adkins
    sadkins@som.umaryland.edu

=cut

use strict;
use warnings;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Pod::Usage;
use LGT::LGTFinder;

############# GLOBALS AND CONSTANTS ################
my $debug = 1;
my ($ERROR, $WARN, $DEBUG) = (1,2,3);
my $logfh;

my $MAX_OVERLAP_THRESHOLD = 20;
my $MIN_LENGTH_THRESHOLD = 0;
####################################################

my %options;

my $results = GetOptions (\%options,
                         "input_file|i=s",
						 "input_list|I=s",
                         "output_dir|o=s",
						 'output_prefix|p=s',
						 'max_overlap|M=i',
						 'min_length|m=i',
						 'ref_lineage|r=s',
                         "log|l=s",
                         "debug|d=s",
                         "help|h"
                          );

&check_options(\%options);

$options{max_overlap} = $MAX_OVERLAP_THRESHOLD if (! $options{max_overlap});
$options{min_length} = $MIN_LENGTH_THRESHOLD if (! $options{min_legth});

my $files = LGT::LGTFinder->findLGT({
		'input_file_list'	=> $options{input_list},
		'output_prefix' => $options{output_prefix},
		'output_dir'	=> $options{output_dir},
		'blast'			=> $options{input_file},
		'ref_lineage'		=> $options{ref_lineage},
		'max_overlap'		=> $options{max_overlap},
		'min_length'		=> $options{min_length}
	});

exit(0);

sub check_options {
   my $opts = shift;
   if( $opts->{'help'} ) {
       &_pod;
   }

   if( $opts->{'log'} ) {
       open( $logfh, "> $opts->{'log'}") or die("Can't open log file ($!)");
   }

   $debug = $opts->{'debug'} if( $opts->{'debug'} );

   foreach my $req ( qw( output_dir ref_lineage ) ) {
       &_log($ERROR, "Option $req is required") unless( $opts->{$req} );
   }

   if (! $opts->{input_file} && ! $opts->{input_list}){
		&_log($ERROR, "Must provide either an --input_file or an --input_list");
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
