#!/usr/bin/env perl
=head1 NAME

lgt_create_validated_bam.pl - Validated a create a BAM file after finding LGT hits

=head1 SYNOPSIS

 USAGE: lgt_create_validated_bam.pl
       --input_file=/path/to/some/input.file
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

B<--clone_hits,-c>
	File that gives the donor/host hits and their LCA lineages.  Obtained from lgt_finder.pl

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

 	This is a wrapper script that will call LGT::LGTSeek->validate_bam to create a BAM file after finding LGT in lgt_finder.pl

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

my $results = GetOptions (\%options,
					 "input_file|i=s",
					 "clone_hits|c=s",
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
	
# Returns a hash of group assignemt counts and the files, but not really necessary to bring out here.
my $val_data = $lgt_obj->validated_bam( {
		'input'		=> $options{input_file},
		'output_dir'=> $options{output_dir},
		'by_clone'	=> $options{clone_hits},
	} );

## Add numbers for validated-LGT to post_processing.tab
my (@header, @vals);
push( @header, 'lgt_valid_blast' );
push( @vals,   "$val_data->{count}" );
LGT::LGTSeek->print_tab( "$output_dir/$name\_validations.tab", \@header, \@vals );

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

   foreach my $req ( qw(input_file clone_hits samtools_path output_dir) ) {
       &_log($ERROR, "Option $req is required") unless( $opts->{$req} );
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
