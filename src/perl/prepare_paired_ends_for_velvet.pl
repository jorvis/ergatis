#!/usr/bin/env perl

=head1 NAME

prepare_input_for_velvet.pl - Will read in a file and run shuffle sequences (provided by velvlet)
    for any paired end. Will simply copy non-paired end data into the output directory

=head1 SYNOPSIS

 USAGE: prepare_short_paired_input_for_velvet.pl
       --output_directory=/path/to/output_dir/
       --short_input_list=/path/to/paired_end_files.list
       --short_format=[fasta|fastq]
       --short_output_list=/path/to/short_reads.list
       --long_input_list=/path/to/long_paired_ends.list
       --long_format=[fasta|fastq]
       --long_output_list=/path/to/long_reads.list
       --velvet_path=/path/to/velvet_dir/
     [ --log=/path/to/file.log
       --debug=3
       --help
     ]

=head1 OPTIONS

B<--output_directory,-o>
    Path to write the shuffled sequence files to.

B<--short_input_list,-s>
    Value should be a list (or comma separated list of lists) that contains 2 files (in the same format)
    which represent paired end files. Accepts only fasta or fastq format.

B<--short_format,-s>
    Specify the format that each of the input lists is. For example if you specified this for -i:
      /path/to/fastq.list,/path/to/fasta.list

    You would specificy -f fastq,fasta for this option.

B<--short_output_list,-so>
    Will create a list of the short read output files

B<--long_input_list,-l>
    Same as short input list, but these are for long reads (ex sanger, 454, etc.)

B<--long_format,-lf>
    Same as short format but for the long input

B<--long_output_list,-lo>
    Sames as for short output list, but with long reads

B<--velvet_path,-v>
    Path to velvet install dir

B<--log,-L>
    Logfile.

B<--debug,-d>
    1,2 or 3. Higher values more verbose.

B<--help,-h>
    Print this message

=head1  DESCRIPTION

    The only reason the short and long need to be separated so we can make two output lists of 
    shuffled paired files so that these can be passed into run_velvet_optimiser.pl script separately.
 
=head1  INPUT
    Input should be a list of two files. The two files should represent one runs paired end data. This script
    will run shuffleSequences.pl on it so it can be used as input into velvet.

=head1 OUTPUT
    Will write two output lists, long_reads.list and short_reads.list.

=head1  CONTACT

    Kevin Galens
    kgalens@gmail.com

=cut

use strict;
use warnings;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Pod::Usage;
use Data::Dumper;

############# GLOBALS AND CONSTANTS ################
my $debug = 1;
my ($ERROR, $WARN, $DEBUG) = (1,2,3);
my $logfh;
my @long_lists;
my @short_lists;
my $velvet_path;
my $out_dir;
####################################################

my %options;
my $results = GetOptions (\%options,
                          "output_directory|o=s",
                          "short_input_list|s=s",
                          "short_format|sf=s",
                          "short_output_list|so=s",
                          "long_input_list|l=s",
                          "long_format|lf=s",
                          "long_output_list|lo=s",
                          "velvet_path|v=s",
                          "log|l=s",
                          "debug|d=s",
                          "help|h"
                          );

&check_options(\%options);

my $fasta_exe = $velvet_path."/shuffleSequences_fasta.pl";
my $fastq_exe = $velvet_path."/shuffleSequences_fastq.pl";

#create the output lists
my $long_output_list = $options{'long_output_list'} || $out_dir."/long_reads.list";
my $lofh;
open($lofh, "> $long_output_list") or die("Can't write to $long_output_list: $!");
my $short_output_list = $options{'short_output_list'} || $out_dir."/short_reads.list";
my $sofh;
open($sofh, "> $short_output_list") or die("Can't write to $short_output_list: $!");

my $bases = {}; #incase we have multiple files of the same name
&shuffle_sequences( $sofh, @short_lists );
&shuffle_sequences( $lofh, @long_lists );

close($sofh);
close($lofh);

sub shuffle_sequences {
    my ($fh, @lists) = @_;
    
    if( !defined( $fh ) ) {
        die("Output list file handle was not passed into sub");
    }


    foreach my $l ( @lists ) {
        my $base;
        if( $l->{'files'}->[0] =~ m|/([^/]+)\.[^\.]+$| ) {
            $base = $1;
        } else {
            die("Couldnt parse basename from $l->{'files'}->[0]");
        }

        if( exists( $bases->{$base} ) ) {
            $bases->{$base}++;
            $base .= ".$bases->{$base}";
        }

        my $out_file = $out_dir."/$base.$l->{'format'}";

        my $exe;
        if( $l->{'format'} eq 'fasta' ) {
            $exe = $fasta_exe;
        } elsif( $l->{'format'} eq 'fastq' ) {
            $exe = $fastq_exe;
        } else {
            die("Bad format: $l->{'format'}");
        }

        my $file_string = join(" ", @{$l->{'files'}} );
        $exe .= " $file_string $out_file";

        system($exe);

        print $fh "$out_file\n";

    }

    
}


sub check_options {

   my $opts = shift;

   if( $opts->{'help'} ) {
       &_pod;
   }

   if( $opts->{'log'} ) {
       open( $logfh, "> $opts->{'log'}") or die("Can't open log file ($!)");
   }

   foreach my $req ( qw(output_directory velvet_path) ) {
       &_log($ERROR, "Option $req is required") unless( $opts->{$req} );
   }

   $out_dir = $opts->{'output_directory'};
   $velvet_path = $opts->{'velvet_path'};
   
   if( $opts->{'short_input_list'} ) {
       die("Option --short_format is required with --short_input_list")
           unless( $opts->{'short_format'} );
       
       @short_lists = &store_files( $opts->{'short_input_list'}, $opts->{'short_format'} );
   }

   if( $opts->{'long_input_list'} ) {
       die("Option --long_format is required with --long_input_list")
           unless( $opts->{'long_format'} );
       
       @long_lists = &store_files( $opts->{'long_input_list'}, $opts->{'long_format'} );
   }
}

sub store_files {
    my ($lists, $formats) = @_;

    my @f = split(/[,\s]+/, $formats );
    my @ls = split(/[,\s]+/, $lists );
    die("The number of items in --short_input_list and --short_format should be equal")
        unless( scalar(@ls) == scalar(@f) );

    my @retval;
    
    for( my $i = 0; $i < @ls; $i++ ) {
        my $list_file = $ls[$i];
        open(LIST, "< $list_file") or die("Couldn't open file $list_file");
        chomp( my @tmp = <LIST> );
        close(LIST);
        die("There should be 2 files in list: $list_file. Found ".scalar(@tmp) )
            unless( @tmp == 2 );
        push(@retval, { "files" => \@tmp, "format" => $f[$i] } );
    }
    return @retval;
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
