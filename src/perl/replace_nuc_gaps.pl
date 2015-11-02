#!/usr/bin/env perl

=head1 NAME

perl_template.pl - Description

=head1 SYNOPSIS

 USAGE: perl_template.pl
       --input_file=/path/to/some/input.file
       --output=/path/to/transterm.file
       --nucleotide=T
     [ --log=/path/to/file.log
       --debug=3
       --help
     ]

=head1 OPTIONS

B<--input_file,-i>

B<--output_file,-o>

B<--nucleotide,-n>
    T or F.  Indicates if FASTA file contains nucleotide sequences (T) or polypeptide sequences (F)

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

my $nuc_flag;
my $in_file;
my $out_file;
my $header;
####################################################

my %options;
my $results = GetOptions (\%options,
                         "input_file|i=s",
                         "output_file|o=s",
                         "nucleotide|n=s",
                         "log|l=s",
                         "debug|d=s",
                         "help|h"
                          );

&check_options(\%options);

if ($nuc_flag) {
	open FILE, $in_file or _log($ERROR, "Cannot open input file $in_file for reading");
	open OUT, ">" . $out_file or _log($ERROR, "Cannot open output file $out_file for writing");
	while (<FILE>){
		my $line = $_;
		next if ($line =~ /^\s+$/);
		if ($line =~ /^>(.+)/) {
			$header = $1;	# Get header name for logging
		} else {
			if ($line =~ /-/) {	# if a dash is encountered, replace with N
				_log($DEBUG, "Replacing '-' with 'N' in sequence for $header");
				$line =~ s/-/N/g
			}
		}
		print OUT $line;
	}
	close OUT;
	close FILE;
# open file 
# read line
# if not >
	# s/-/N/g
	# write to output

} else {
	_log($DEBUG, "File contains polypeptide sequences.  Copying input file to output directory");
	system("cp $in_file out_file")
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

   foreach my $req ( qw(input_file output_file) ){
       &_log($ERROR, "Option $req is required") unless( $opts->{$req} );
   }
   
   if (! exists($opts->{'nucleotide'}) ){
   	$nuc_flag =0;
   	&_log($WARN, "Wasn't clear if FASTA sequences were nucleotide or polypeptide.  Will not change gaps.");
   }	#Take the safe route and not change anything if nothing was specified
   $nuc_flag = (uc($opts->{'nucleotide'}) eq "T") ? 1 : 0;
   $in_file = $opts->{'input_file'};
   $out_file = $opts->{'output_file'};
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
