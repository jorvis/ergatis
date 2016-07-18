#!/usr/bin/env perl
=head1 NAME

lgt_mpileup.pl - Wrapper to call samtools 'mpileup' method

=head1 SYNOPSIS

 USAGE: lgt_mpileup.pl
       --input_file=/path/to/some/input.file
	   --fasta_ref=/path/to/reference.fasta
       --output_dir=/path/to/transterm
	   --clone_hits=/path/to/hits_by_clone.txt
	   --samtools_path=/usr/local/bin/samtools
     [ --log=/path/to/file.log
       --debug=3
       --help
     ]

=head1 OPTIONS

B<--input_file,-i>
	A BAM file that has been filtered by Picard-tools and Prinseq

B<--fasta_ref,-f>
	A nucleotide fasta reference file.  Can also accept a list of references (with a .list extension).

B<--output_dir,-o>

B<--samtools_path,-s>
    Path to samtools executable

B<--log,-l>
    Logfile.

B<--debug,-d>
    1,2 or 3. Higher values more verbose.

B<--help,-h>
    Print this message

=head1  DESCRIPTION

 	This is a wrapper script that will call LGT::LGTSeek->mpileup to create a .pileup format file for finding potential SNPs

=head1  INPUT

=head1 OUTPUT

=head1  CONTACT

    Shaun Adkins
    sadkins@som.umaryland.edu

=cut

use strict;
use warnings;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Pod::Usage;
use File::Basename;
use LGT::LGTSeek;

############# GLOBALS AND CONSTANTS ################
my $debug = 1;
my ($ERROR, $WARN, $DEBUG) = (1,2,3);
my $logfh;
####################################################
my $output_dir;
my %options;
my $count_id = 0;
my @ref_files;

my $results = GetOptions (\%options,
					 "input_file|i=s",
					 "fasta_ref|f=s",
                     "samtools_path|s=s",
                     "output_dir|o=s",
                     "log|l=s",
                     "debug|d=s",
                     "help|h"
                      );

&check_options(\%options);

my $lgt_obj = LGT::LGTSeek->new( {
		'samtools_bin' => $options{samtools_path},
		'verbose' => 1
} );

foreach my $ref (@ref_files){
	chomp $ref;
	my $mpilup_output = $lgt_obj->mpileup( {
		    'input'      => $options{input_file},
	        'output_dir' => $output_dir,
	        'ref'        => $ref
	    }
	);
}

exit(0);

# Process read options
sub check_options {
   my $opts = shift;
   if( $opts->{'help'} ) {
       &_pod;
   }

   if( $opts->{'log'} ) {
       open( $logfh, "> $opts->{'log'}") or die("Can't open log file ($!)");
   }

   $debug = $opts->{'debug'} if( $opts->{'debug'} );

   foreach my $req ( qw(input_file fasta_ref samtools_path output_dir) ) {
       &_log($ERROR, "Option $req is required") unless( $opts->{$req} );
   }


   # Gather reference files together
   if($opts->{fasta_ref} =~ /\.list$/) {
       @ref_files = `cat $opts->{fasta_ref}`;
   } else {
       push @ref_files, $opts->{fasta_ref};
   }

   $output_dir = $opts->{'output_dir'};
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
