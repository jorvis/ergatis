#!/usr/bin/perl

=head1 NAME

run_velvet_optimiser.pl - Eragtis wrapper script for running Velvet Optimiser.

=head1 SYNOPSIS

 USAGE: run_velvet_optimiser.pl
       --input_file=/path/to/file.reads
       --velvet_path=/abs/path/to/velvet_dir
       --file_format=-fasta
       --read_type=-short
       --start_hash_length=17
       --end_hash_length=31
       --other_optimiser_opts='-a yes'
       --output_directory=/path/to/out
       --output_list=/path/to/contig.fa.list
     [ --log=/path/to/file.log
       --debug=3
       --help
     ]

=head1 OPTIONS

B<--input_file,-i>
    Path to an input file to be passed into the optimiser

B<--velvet_path,-p>
    Absolute directory of where the velvet executables are held

B<--file_format,-f>
    Can be one of the following:
    [-fasta|-fastq|-fasta.gz|-fastq.gz|-eland|-gerald]

B<--read_type,-r>
    Can be one of the following:
    [-short|-shortPaired|-short2|-shortPaired2|-long|-longPaired]

B<--start_hash_length,-s>
    The hash length (or kmer) for the first iteration of velvet.
    Must be an odd number

B<--end_hash_length,-e>
    The hash length (or kmer) for the last iteration of velvet
    Must be an odd number

B<--other_optimiser_opts,-o>
    Right now, should only be -a yes.  But will be passed directly
    into the executable. So don't break things and such.

B<--output_directory,-O>
    Directory to write output

B<--output_list,-u>
    The path to the list file which will contain the full path
    to the contigs.fa produced by velvet.  Optional; if not provided,
    no list will be made.

B<--log,-l>
    Logfile.

B<--debug,-d>
    1,2 or 3. Higher values more verbose.

B<--help,-h>
    Print this message

=head1  DESCRIPTION

    Will take in an input file and run velvet optimiser.
 
=head1  INPUT
    Anything valid that velvet can handle.

=head1 OUTPUT

=head1  CONTACT

    Kevin Galens
    kgalens@gmail.com

=cut

use strict;
use warnings;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Pod::Usage;
use File::Find;
$|++;

############# GLOBALS AND CONSTANTS ################
my $debug = 1;
my ($ERROR, $WARN, $DEBUG) = (1,2,3);
my $logfh;
my @input_files;
####################################################

my %options;
my $results = GetOptions (\%options,
                          "input_file|i=s",
                          "velvet_path|v=s",
                          "file_format|f=s",
                          "read_type|r=s",
                          "start_hash_length|s=s",
                          "end_hash_length|e=s",
                          "other_optimiser_opts|o=s",
                          "output_directory|O=s",
                          "output_list|u=s",
                          "log|l=s",
                          "debug|d=s",
                          "help|h"
                          );

&check_options(\%options);

#find the necessary executables
my $velvet_optimiser_exec = $options{'velvet_path'}."/contrib/VelvetOptimiser/VelvetOptimiser.pl";
&_log($ERROR, "Could not find VelvetOptimiser.pl (executable). Thought it would be here: ".
      "$velvet_optimiser_exec.") unless( -e $velvet_optimiser_exec );

#create the command string
my $input_file_string = join(" ", @input_files);
my $opts_string = "$options{'read_type'} $options{'file_format'} $input_file_string";
my $cmd = "$velvet_optimiser_exec -f '$opts_string' -s $options{'start_hash_length'} -e ".
    "$options{'end_hash_length'}";
$cmd .= " ".$options{'other_optimiser_opts'} if( $options{'other_optimiser_opts'} );
$cmd .= " 2>\&1";

#Before running the command, switch to the correct output directory
#since it will just write the output to that directory.
$cmd = "cd $options{'output_directory'}; $cmd";

#Also need to add velvet_dir to the path
$cmd = "export PATH=\$PATH:$options{'velvet_path'}; $cmd";

#And run it, making sure to pass along stdout
&_log($DEBUG, $cmd);
my $flag = 0;
my $data_dir;
open(CMD, "$cmd |") or &_log($ERROR, "Could not run command $cmd ($!)");
while( <CMD> ) {
    chomp;

    #This indicates we found the indicator line previously and
    #this line should contain the directory to where the output is
    if( $flag ) {
        $data_dir = $_;
        if( !-d $data_dir ) {
            &_log($ERROR, "Parsed data_dir: $data_dir and it does not exist");
        }
        $flag = 0;

    #When we come across this line, it means the output directory
    #is next.
    } elsif( $_ =~ /Assembly output files are in the/ ) {
        $flag = 1;
    }

    #We print all lines that VelvetOptimiser produces to the log
    #file and also to STDOUT (only if debug level high enough).
    &_log($DEBUG, $_);
}
close(CMD);

#make the output list if the user specified a place to make it.
if( $options{'output_list'} ) {
    #we should make a list file
    open( my $out, "> $options{'output_list'}" ) or die("Can't open $options{'output_list'} $!");
    
    #Make sure it exists. If not die.
    my $contig_file = $data_dir."/contigs.fa";
    if( !-e $contig_file ) {
        &_log($ERROR, "Could not find a contigs.fa file in the data directory $data_dir");
    }
    print $out $contig_file."\n";
    close($out);
}


sub check_options {

   my $opts = shift;

   if( $opts->{'help'} ) {
       &_pod;
   }

   mkdir("$opts->{'output_directory'}") unless( -d $opts->{'output_directory'} );

   if( $opts->{'log'} ) {
       open( $logfh, "> $opts->{'log'}") or die("Can't open log file ($opts->{'log'}) ($!)");
   }

   foreach my $req ( qw(input_file output_directory velvet_path file_format 
                        read_type start_hash_length end_hash_length) ) {
       &_log($ERROR, "Option $req is required") unless( $opts->{$req} );
   }

   open( TEST, "< $opts->{'input_file'}" ) or &_log( $ERROR, "Could not open $opts->{'input_file'} ($!)");
   chomp( my @lines = <TEST> );
   my $first_line = $lines[0];
   if( $first_line =~ m|^/| ) {
       #must be a list file
       @input_files = @lines;
   } else {
       push(@input_files, $opts->{'input_file'});
   }
   close(TEST);

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
