#!/usr/bin/perl

use strict;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Pod::Usage;
use POSIX;

my %options = ();
my $results = GetOptions (\%options,
                          'first_part_dir|f=s',
			  'second_part_dir|s=s', 
                          'out_dir|o=s',
			  'log|l=s',
                          'debug=s',
                          'help|h') || pod2usage();

## display documentation
if( $options{'help'} ){
    pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDERR} );
}

$| = 1;

my $firstPartitionSetDir = $options{first_part_dir};

opendir DOT, $firstPartitionSetDir or die("Can't open $firstPartitionSetDir for scanning.\n");

my @firstPartFiles = sort { $a cmp $b } map { "$firstPartitionSetDir/$_" } grep { /^yy04_partition_1_groups.group\d+\.fa$/ } readdir DOT;

closedir DOT;

my $secondPartitionSetDir = $options{second_part_dir};

opendir DOT, $secondPartitionSetDir or die("Can't open $secondPartitionSetDir for scanning.\n");

my @secondPartFiles = sort { $a cmp $b } map { "$secondPartitionSetDir/$_" } grep { /^yy07_part2_filtered.group\d+\.fa$/ } readdir DOT;

closedir DOT;

my @finalFileSet = ( @firstPartFiles, @secondPartFiles );

my $counter = 0;

my $outputDir = $options{out_dir};

foreach my $file ( @finalFileSet ) {
   
   my $printCounter = &printable($counter);
   
   my $target = "$outputDir/input_$printCounter.fa";

   system("ln -sf $file $target");

   if ( not -l $target ) {
      
      print "input_linker.pl failed to create link \"$target\".  Aborting.\n";

      exit(1);
   }

   $counter++;
}

print "input_linker.pl completed successfully.\n";

exit(0);

sub printable {
   
   my $num = shift;

   if ( $num < 10 ) {
      
      return '00' . $num;

   } elsif ( $num < 100 ) {
      
      return '0' . $num;

   } else {
      
      return $num;
   }
}
