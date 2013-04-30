#!/usr/bin/perl

use strict;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Pod::Usage;
use POSIX;

my %options = ();
my $results = GetOptions (\%options,
			  'input_dir|i=s', 
                          'log|l=s',
                          'debug=s',
                          'help|h') || pod2usage();


## display documentation
if( $options{'help'} ){
    pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDERR} );
}

#$| = 1;

my $pwd = $options{input_dir};

opendir DOT, "$pwd" or die("Can't open $pwd for scanning.\n");

my @files = map { "$pwd/$_" } grep { /^yy04_partition_1_groups\.group\d+\.fa$/ } readdir DOT;

closedir DOT;

my $finalFile = '';

foreach my $file ( sort { $a cmp $b } @files ) {
   
   $finalFile = $file;
}

my $targetFile = $finalFile;


$targetFile =~ s/\/([^\/]+)$/\//;

$targetFile .= 'yy05_partition_1_groups.finalGroup.fa';
#print "targetFile: $targetFile\nfinalFile:  $finalFile\n";

system("rm -f $targetFile");
system("mv $finalFile $targetFile");

if ( -e $targetFile ) {
   
   print "rename_final_group.pl completed successfully.\n";

   exit(0);

} else {
   
   print "rename_final_group.pl failed.\n";

   exit(1);
}

