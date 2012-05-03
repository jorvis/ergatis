#!/usr/bin/env perl

=head1 NAME

mpileup_wrapper.pl - Description

=head1 SYNOPSIS

 USAGE: mpileup_wrapper.pl
       --tmp_directory=/tmp/directory
       --mpileup_options="-C50 -b whatever.bam"
       --fasta_file=/path/to/fasta

=head1  DESCRIPTION

  Since the fasta file needs to be indexed and the index needs to be located where
  the fasta file is located, this script will make symlinks to the fasta files in
  the temp directory and indexes will be created there. This solves an issue about
  using input fasta files when you don't have write access to the directory (or having
  index files around in unrelated directories)

  Any other options should be specified in the mpileup options input parameter (including
  the input bam files). The bam files are expected to be sorted/indexed.
 
=head1  CONTACT

    Kevin Galens
    kgalens@gmail.com

=cut

use strict;
use warnings;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use File::Basename;

############# GLOBALS AND CONSTANTS ################
my $debug = 1;
my ($ERROR, $WARN, $DEBUG) = (1,2,3);
my $logfh;
####################################################

my %options;
my $results = GetOptions (\%options,
						  "tmp_directory|t=s",
						  "fasta_file|f=s",
						  "mpileup_options|m=s",
						  "samtools_exec|s=s",
                          );

&check_options(\%options);

## Create symlink
my $basename = basename( $options{'fasta_file'} );
my $new_file = $options{'tmp_directory'}."/$basename";
symlink( $options{'fasta_file'}, $new_file );

## mpileup will automatically create index if it doesn't exist.
my $cmd = $options{'samtools_exec'}." mpileup $options{'mpileup_options'} -f $new_file";
system($cmd) && die("Problem running command $cmd. Exit val: $?");

sub check_options {
   my $opts = shift;

   foreach my $req ( qw(tmp_directory fasta_file mpileup_options samtools_exec) ) {
       &_log($ERROR, "Option $req is required") unless( $opts->{$req} );
   }
}
