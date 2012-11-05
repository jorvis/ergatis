#!/usr/bin/env perl

=head1 NAME

gkpstore2fastq.pl - Description

=head1 SYNOPSIS

 USAGE: gkpstore2fastq.pl
       --input_gkpstore=/path/to/FILE.gkpStore
       --output_fastq=/path/to/output.fastq
     [ --log=/path/to/file.log
       --debug=3
       --help
     ]

=head1  DESCRIPTION

 DESCRIPTION

=head1  CONTACT

    Kevin Galens
    kgalens@gmail.com

=cut

use strict;
use warnings;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Pod::Usage;

############# GLOBALS AND CONSTANTS ################
my $debug = 1;
my ($ERROR, $WARN, $DEBUG) = (1,2,3);
my $logfh;
my %gexes; # stores various versions of the gatekeeper executable
####################################################

my %options;
my $results = GetOptions (\%options,
			  "input_gkpstore|i=s",
			  "output_directory|o=s",
			  "gatekeeper_exec_6-1=s",
			  "gatekeeper_exec_7-0=s",
			  "output_prefix|p=s",
			  "log|l=s",
			  "debug|d=s",
			  "help|h"
                          );

&check_options(\%options);

my $gkpstore_version = $1 if( $options{'input_gkpstore'} =~ /_CA([\d\.]+)/ );
die("Could not find gatekeeper exec version [$gkpstore_version] for $options{'input_gkpstore'}")
    unless( exists( $gexes{$gkpstore_version} ) );

my $gexe = $gexes{$gkpstore_version};
my $cmd = "$gexe -dumpfastq $options{'output_directory'}/$options{'output_prefix'} $options{'input_gkpstore'}";
system($cmd) && die("Could not run command [$cmd]");

sub check_options {

   my $opts = shift;

   if( $opts->{'help'} ) {
       &_pod;
   }

   if( $opts->{'log'} ) {
       open( $logfh, "> $opts->{'log'}") or die("Can't open log file ($!)");
   }

   $debug = $opts->{'debug'} if( $opts->{'debug'} );

   foreach my $req ( qw(input_gkpstore output_directory output_prefix) ) {
       &_log($ERROR, "Option $req is required") unless( $opts->{$req} );
   }
   

   foreach my $gexe ( qw(gatekeeper_exec_6-1 gatekeeper_exec_7-0) ) {
       next unless( exists( $opts->{$gexe} ) );
       my $version = $1 if( $gexe =~ /\_([\d\-]+)$/ );
       $version =~ s/-/\./g;
       $gexes{$version} = $opts->{$gexe};
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
