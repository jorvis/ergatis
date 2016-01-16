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

B<--hits_list,-c>
	List file that lists the clone and trace files containing the donor/host hits and their LCA lineages.
	Should contain just one trace and clone hit file each. Obtained from lgt_finder.pl

B<--output_dir,-o>

B<--prefix, -p>
	Name you would like to give as a prefix in the output file name

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
use File::OpenFile qw(open_file);
use LGT::LGTSeek;

############# GLOBALS AND CONSTANTS ################
my $debug = 1;
my ($ERROR, $WARN, $DEBUG) = (1,2,3);
my $logfh;

my $PREFIX = "validation";
####################################################
my $output_dir;
my %options;
my $count_id = 0;
my $name;

my $results = GetOptions (\%options,
					 "input_file|i=s",
					 "hits_list|h=s",
                     "samtools_path|s=s",
                     "output_dir|o=s",
					 "prefix|p=s",
                     "log|l=s",
                     "debug|d=s",
                     "help|h"
                      );

&check_options(\%options);

my $hits = parse_list($options{hits_list});

my $lgt_obj = LGT::LGTSeek->new( {
		'samtools_bin' => $options{samtools_path},
		'verbose' => 1
} );

# Returns a hash of group assignemt counts and the files, but not really necessary to bring out here.
my $val_data = $lgt_obj->validated_bam( {
		'input'		=> $options{input_file},
		'output_dir'=> $options{output_dir},
		'output_prefix' => $name,
		'by_clone'	=> $hits->{by_clone}
	} );

## Add numbers for validated-LGT to post_processing.tab
my (@header, @vals);
push( @header, 'lgt_valid_blast' );
push( @vals,   "$val_data->{count}" );
LGT::LGTSeek->print_tab( "$output_dir/$name\_validations.tab", \@header, \@vals );

exit(0);

# Parse the list file to get the relevant lists.
sub parse_list {
	my $input_list = shift;
	my %hits;
	&_log($DEBUG, "Parsing $input_list\n");
    my $oh = open_file( $input_list, "in" );
	while (<$oh>) {
		chomp;
		if (/by_clone/){
			&_log($ERROR, "Multiple clone hits file found in input list file") if ($hits{by_clone});
			$hits{by_clone} = $_;
		}
		elsif (/by_trace/){
			&_log($ERROR, "Multiple trace hits file found in input list file") if ($hits{by_trace});
			$hits{by_trace} = $_;
		}
	}

	return \%hits;
}


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

   foreach my $req ( qw(input_file hits_list samtools_path output_dir) ) {
       &_log($ERROR, "Option $req is required") unless( $opts->{$req} );
   }

   $output_dir = $opts->{'output_dir'};

   $name = $opts->{prefix} ? $opts->{prefix} : $PREFIX;
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
