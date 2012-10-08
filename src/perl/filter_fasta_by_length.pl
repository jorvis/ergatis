#!/usr/bin/perl

eval 'exec /usr/bin/perl  -S $0 ${1+"$@"}'
    if 0; # not running under some shell

eval 'exec /usr/bin/perl  -S $0 ${1+"$@"}'
    if 0; # not running under some shell
BEGIN{foreach (@INC) {s/\/usr\/local\/packages/\/local\/platform/}};
use lib (@INC,$ENV{"PERL_MOD_DIR"});
no lib "$ENV{PERL_MOD_DIR}/i686-linux";
no lib ".";


=head1 NAME

filter_fasta_by_length.pl - Remove sequences shorter than specified length
                             
=head1 SYNOPSIS

USAGE: ./filter_fasta_by_length.pl --limit=##
				   --input=/path/to/input/file
				   --output=/path/to/output/file

=cut

use strict;
use warnings;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Pod::Usage;
use File::Basename;
use Ergatis::Logger;


$| = 1;

my $logger;
my %options = &_parse_options();
my $input_file = $options{'input_file'};
my $output_file = $options{'output_file'};
my $lengthCutoff = $options{'limit'};


open (IN, "<$input_file") or $logger->logdie("Can't open $input_file for reading.\n");
open (OUT, ">$output_file") or $logger->logdie("Could not write to $output_file: $!");

my $first = 1;
my $currentHeader = '';
my $currentSeq = '';


while ( my $line = <IN> ) {
   
   if ( $line =~ /^>/ ) {
      
      if ( $first ) {
	 
	 $currentHeader = $line;

	 $first = 0;

      } else {
	 
	 my $currentLength = length($currentSeq);

	 if ( $currentLength >= $lengthCutoff ) {
	    
	    print OUT $currentHeader;

	    print OUT &fastaFy($currentSeq);
	 }

         $currentHeader = $line;

         $currentSeq = '';
      }

   } else {
      
      chomp $line;

      $currentSeq .= $line;
   }
}

close OUT;
close IN;

exit(0);


###subs

sub fastaFy {

   my $seq = shift;

   my $index = 0;

   my $result = '';

   my $lineLength = 60;

   while ( length($seq) - $index > $lineLength ) {

      $result .= substr($seq, $index, $lineLength) . "\n";

      $index += $lineLength;
   }

   $result .= substr($seq, $index) . "\n";

   return $result;
}


# Parse command-line arguments                                       
sub _parse_options {
    my %opts = ();

    GetOptions(\%opts,
                'input_file|i=s',
                'output_file|o=s',
                'limit|l=i',
                'help' ) || pod2usage();

    if ($opts{'help'}) {
        pod2usage ( { -exitval => 0, -verbose => 2, -output => \*STDERR } );
    }
    
    my $logfile = Ergatis::Logger::get_default_logfilename();
    my $debug = 4;
    $logger = new Ergatis::Logger( 'LOG_FILE'   =>  $logfile,
                                   'LOG_LEVEL'  =>  $debug );
    $logger = Ergatis::Logger::get_logger();
    

    defined ($opts{'output_file'}) || $logger->logdie("Please specify an output file");
    defined ($opts{'input_file'}) || $logger->logdie("Please specify an input file");
    defined ($opts{'limit'}) || $logger->logdie("Please specify a limit");
    

    return %opts

}




