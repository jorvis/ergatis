#!/usr/local/bin/perl

=head1 NAME

create_clovrtag_iterator_list.pl - Description

=head1 SYNOPSIS

 USAGE: create_clovrtag_iterator_list.pl
       --input_tags=TAG1,TAG2,TAG3
       --output_iter_list=/path/to/output_iter_list
     [ --log=/path/to/file.log
       --debug=3
       --help
     ]

=head1 OPTIONS


=head1  DESCRIPTION

    Need to make these vars: $;INPUT_TAG_NAME$;
 
=head1  INPUT
    Describe the input

=head1 OUTPUT
    Describe the output

=head1  CONTACT

    Kevin Galens
    kgalens@gmail.com

=cut

use strict;
use warnings;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Pod::Usage;
use File::OpenFile qw(open_file);

############# GLOBALS AND CONSTANTS ################
my $debug = 1;
my ($ERROR, $WARN, $DEBUG) = (1,2,3);
my $logfh;
my $outlist;
my @tags;
####################################################

my %options;
my $results = GetOptions (\%options,
                          "input_tags|i=s",
                          "output_iter_list|o=s",
                          "log|l=s",
                          "debug|d=s",
                          "help|h"
                          );

&check_options(\%options);

my $ofh = open_file( $outlist, "out");
print $ofh '$;INPUT_TAG_NAME$;'."\n";
map { print $ofh "$_\n"; } @tags;
close($ofh);

sub check_options {

   my $opts = shift;

   if( $opts->{'help'} ) {
       &_pod;
   }

   if( $opts->{'log'} ) {
       open( $logfh, "> $opts->{'log'}") or die("Can't open log file ($!)");
   }

   foreach my $req ( qw(input_tags output_iter_list) ) {
       &_log($ERROR, "Option $req is required") unless( $opts->{$req} );
   }

   $outlist = $opts->{'output_iter_list'};
   @tags = split(/[,\s]+/, $opts->{'input_tags'} );
   die("Could not parse any tags from input") unless( @tags > 0 );
}

sub _log {
   my ($level, $msg) = @_;
   if( $level <= $DEBUG ) {
      print STDOUT "$msg\n";
      print $logfh "$msg\n" if( defined( $logfh ) );
      exit(1) if( $level == $ERROR );
   }
}

sub _pod {
    pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDERR} );
}
