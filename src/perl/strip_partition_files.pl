#!/usr/bin/perl

use strict;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Pod::Usage;
use POSIX;
use File::Basename;
#BEGIN {
#use Ergatis::Logger;
#}

my %options = ();
my $results = GetOptions (\%options,
                          'input_dir|i=s',
			  'output_dir|o=s',
			  'exec_dir|e=s', 
			  'packages_dir|p=s',
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

my $pwd = $options{input_dir};
my $out_dir = $options{output_dir};
my $exec = $options{exec_dir};
my $packages_dir = $options{packages_dir};

my $linkDir = "$pwd";

opendir DOT, $linkDir or die("Can't open $linkDir for scanning.\n");

my @files = sort { $a cmp $b } map { "$linkDir/$_" } grep { /^input_\d+\.fa$/ } readdir DOT;

closedir DOT;

foreach my $file ( @files ) {
   my $file_base = basename($file);   
   my $command = "PYTHONPATH=$packages_dir/khmer/python $packages_dir/Python-2.7/bin/python2.7 $exec/strip-partition.py $file > $out_dir/$file_base.strip";

   system($command);

   if ( not -e "$out_dir/$file_base.strip" ) {
      
      print "strip_partition_files.pl failed to create file \"$out_dir/$file_base.strip\".  Aborting.\n";

      exit(1);
   }
}

print "strip_partition_files.pl completed successfully.\n";

