#!/usr/bin/env perl

use strict;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Pod::Usage;
use POSIX;
#BEGIN {
#use Ergatis::Logger;
#}

my %options = ();
my $results = GetOptions (\%options,
                          'input_dir|i=s', 
         		  'output_dir|o=s',
                          'log|l=s',
                          'debug=s',
                          'help|h') || pod2usage();

#my $logfile = $options{'log'} || Ergatis::Logger::get_default_logfilename();
#my $logger = new Ergatis::Logger('LOG_FILE'=>$logfile,
#                                  'LOG_LEVEL'=>$options{'debug'});
#$logger = $logger->get_logger();

## display documentation
if( $options{'help'} ){
    pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDERR} );
}

$| = 1;

my $inDir = $options{input_dir};
my $outDir = $options{output_dir};

opendir DOT, $inDir or die("Can't open $inDir for scanning.\n");

my @files = map { "$inDir/$_" } sort { $a cmp $b } grep { /strip$/ } readdir DOT;

closedir DOT;

foreach my $file ( @files ) {
   
   $file =~ /input_(\d+)\.fa\.strip/;

   my $indexString = $1;

   my $targetFile = "$outDir/targetIDs_$indexString.txt";

   my $grepCommand = "grep -h \'>\' $file | perl -p -i -e \'s/^>//\' | perl -p -i -e \'s/\\/\\d+\$//\' > $targetFile";

   system($grepCommand);

   if ( not -e $targetFile ) {
      
      print "No $targetFile, get_read_IDs.pl failed.\n";

      exit(1);
   }
}

print "get_read_IDs.pl completed successfully.\n";

exit(0);

