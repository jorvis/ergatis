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
use File::Basename;
use Pod::Usage;

############# GLOBALS AND CONSTANTS ################
my $debug = 1;
my ($ERROR, $WARN, $DEBUG) = (1,2,3);
my $logfh;
my $default_exe = "/local/projects/grc/devel/packages/wgs-cvs/Linux-amd64/bin/gatekeeper";
####################################################

my %options;
my $results = GetOptions (\%options,
			  "input_gkpstore|i=s",
			  "output_directory|o=s",
			  "output_prefix|p=s",
			  "log|l=s",
			  "debug|d=s",
			  "help|h"
                          );

&check_options(\%options);

## This will look for the runCA-logs directory in relation to the gkpstore directory
## If the gkpstore directory has been moved, or the logs deleted, this won't work
## and the default_exe will be used. This probably won't work.
my $gexe = &find_gatekeeper_exe( $options{'input_gkpstore'} );
$gexe = $default_exe unless( $gexe and -e $gexe );

my $cmd = "$gexe -dumpfastq $options{'output_directory'}/$options{'output_prefix'} $options{'input_gkpstore'}";
system($cmd) && die("Could not run command [$cmd]");

sub find_gatekeeper_exe {
    my ($gkpstore) = @_;

    my $dir = dirname( $gkpstore );
    my $calogs_dir = $dir."/runCA-logs";
    $calogs_dir = $dir."/../runCA-logs" unless( -d $calogs_dir );

    return unless( -d $calogs_dir );
    
    opendir( DIR, $calogs_dir ) or die("Can't opendir $calogs_dir: $!");
    my @gatekeeper_logs = map { $_ = $calogs_dir."/$_" } grep { /gatekeepr/ } readdir( DIR );
    close( DIR );

    open( IN, "< $gatekeeper_logs[0]" ) or die("can't open $gatekeeper_logs[0]: $!");
    my $f = 0;
    my $exe;
    while( my $line = <IN> ) {
	if( $line =~ /COMMAND/ ) {
	    $f = 1;
	} elsif( $f && $line =~ /^(\S+)/ ) {
	    $exe = $1;
	}
    }
    close(IN);

    die("Found logfile $gatekeeper_logs[0], but could not parse exe") unless( $exe && -e $exe );

    return $calogs_dir;
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
