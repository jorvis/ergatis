#!/usr/bin/env perl

=head1 NAME

filter_fastq.pl - Will filter input fastq files with input list of read names

=head1 SYNOPSIS

 USAGE: filter_fastq.pl
       --fastq_input=/path/to/fastq
       --read_names=/path/to/read_names.txt
       --output_dir=/path/to/outdir/
     [ --log=/path/to/file.log
       --debug=3
       --help
     ]

=head1 OPTIONS

B<--fastq_input,-i>
    Can be a fastq file or list of fastq files (for paired end).

B<--read_names,-r>
    A file with a list of read_names.

B<--output_dir,-o>
    Path to output directory. Will create an output file for each input fastq file.

B<--log,-l>
    Logfile.

B<--debug,-d>
    1,2 or 3. Higher values more verbose.

B<--help,-h>
    Print this message

=head1  DESCRIPTION

    Developed when an illumina dataset has been sub-sampled before assembly to reduce the 
    coverage (since some assemblers can't handle very deep coverage). We often don't keep
    the actual fastq files used for the assemblies around; instead we just save the read
    names used. This will regenerate the fastq files use for assembly, given the original
    fastq file(s) and the list of read names.

=head1  CONTACT

    Kevin Galens
    kevingalens@gmail.com

=cut

use strict;
use warnings;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Pod::Usage;
use File::Basename;

use Data::Dumper;

############# GLOBALS AND CONSTANTS ################
my $debug = 1;
my ($ERROR, $WARN, $DEBUG) = (1,2,3);
my $logfh;
my @fastqs;
####################################################

my %options;
my $results = GetOptions (\%options,
						  "fastq_input|i=s",
						  "read_names|r=s",
						  "output_dir|o=s",
						  "log|l=s",
						  "debug|d=s",
						  "help|h"
						 );

&check_options(\%options);

my @filtered_fastqs = ();
foreach my $fq ( @fastqs ) {
  my $basename = basename( $fq, qw(.fq .fastq .txt) );
  my $outfile = "$options{'output_dir'}/$basename.fastq";
  &filter_fq_file( $fq, $outfile, $options{'read_names'} );
  push(@filtered_fastqs, $outfile);
}

#Create a list file for the output filtered fastq files.
my $basename = basename( $options{'fastq_input'}, qw(.list .fastq .fq .txt) );
my $outlist = $options{'output_dir'}."/$basename";
open(OUT, "> $outlist");
local $" = "\n";
print OUT "@filtered_fastqs\n";
close(OUT);

sub filter_fq_file {
	my ($old, $new, $seqs_file) = @_;

	my %seqs;
	open(my $sfh, "< $seqs_file") or die("Can't open $seqs_file: $!");
	while( my $line = <$sfh> ) {
	  chomp( $line );

	  # Take everything up to the first / or space. I match some more based on 
	  # assumptions I have about the format of the read IDs which I'm not 100%
	  # sure are true. They aren't really used though.
	  my $read_name = $1 if( $line =~ /^\@?(\S+?)((\/|\s)(1|2)).*?$/ );
	  die("Could not parse readname from $line") unless( $read_name );
	  $seqs{$read_name} = 1;
	}
	close($sfh);

	open(OUT, "> $new") or die("Can't open $new for writing: $!");
	open(IN, "< $old") or die("Can't open $old: $!");

	my $count;

	my $instrument_name;
	my $read_info = [];

	# Setup a subref to process each fastq entry
	my $print_read_info = sub {
		my ($fastq_entry) = @_;

		## parse the id
		my $read_name = $1 if( $fastq_entry->[0] =~ /^\@(\S+?)((\/|\s)(1|2)).*?$/ );
		die("Could not parse read_name from $fastq_entry->[0]") unless( $read_name );
		
		if( exists( $seqs{$read_name} )  ) {
		  local $" = "\n";
		  print OUT "@{$fastq_entry}\n";
		}
		
	};

	while( my $line = <IN> ) {
		chomp( $line );
		
		# Grab instrument name from first line
		if( !defined( $instrument_name ) ) {
			$instrument_name = $1 if( $line =~ /^\@([^\:]+)\:/ );
			die("Could not grab instrument name from first line: $line") unless( $instrument_name );
		} 

		if( @{$read_info} != 0 && $line =~ /^\@$instrument_name/ ) {
			#check to see if there are the correct number of lines
			$print_read_info->($read_info);
			$read_info = [];
		}
		
		push(@{$read_info}, $line);
	}

	$print_read_info->( $read_info );

	close(IN);
	close(OUT);

	return $new;
}

sub check_options {

   my $opts = shift;

   if( $opts->{'help'} ) {
       &_pod;
   }

   if( $opts->{'log'} ) {
       open( $logfh, "> $opts->{'log'}") or die("Can't open log file ($!)");
   }

   foreach my $req ( qw(fastq_input output_dir read_names) ) {
       &_log($ERROR, "Option $req is required") unless( $opts->{$req} );
   }

   # Check to see if the input fastq is a list or files
   open( IN, "< $opts->{'fastq_input'}" ) or die("Could not open $opts->{'fastq_input'}: $!");
   chomp( my $fl = <IN> );
   if( $fl =~ /^\@/ ) {
	 push(@fastqs, $opts->{'fastq_input'} );
   } elsif( -e $fl ) {
	 chomp( my @the_rest = <IN> );
	 push(@fastqs, $fl, @the_rest);
   } else {
	 die("Could not determine format of input file $opts->{'fastq_input'}");
   }
   close(IN);

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
