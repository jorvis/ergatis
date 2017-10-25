#!/usr/bin/env perl

=head1  NAME 

lgt_bwa.pl

=head1 SYNOPSIS


      
=head1 OPTIONS

=over 8

This help message is no good.

=back

=head1   DESCRIPTION


=head1 INPUT



=head1 OUTPUT


=head1 CONTACT

   Kevin Galens
   kgalens@gmail.com

   Adapted from lgt_bwa.pl

=cut
use strict;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use File::Basename;
use Pod::Usage;

my %options = ();
my $results = GetOptions (\%options, 
						  'query_file|q=s',
						  'ref_file|r=s',
						  'output_dir|o=s',
						  'tmp_dir|td=s',

						  'index_bam|ib=s',

						  'bwa_path|b=s',
						  'samtools_path|st=s',

						  'mismatch|mm=s',
						  'max_gaps|mg=s',
						  'max_gap_extensions|mge=s',
						  'open_gap_penalty|og=s',
						  'extend_gap_penalty|eg=s',
						  'threads|t=s',
						  'num_aligns|na=s',

						  'cleanup|c=s',
						  'help|h');

##############################
# Global Vars                #
##############################
my $PAIRED = 0;
my $query_files = [];
my $prefix = "";
my $options_string = '';
##############################

## make sure all passed options are peachy
&check_parameters(\%options);

## are the input reference files indexed?
my $ref_file = $options{'ref_file'};
$ref_file = &create_index( $ref_file ) unless( index_exists( $ref_file ) );

## create bam file
my $bam = &create_bam( $ref_file, $query_files, $prefix );

## sort/index the bam if if necessary
&index_bam( $bam ) if( $options{'index_bam'} );

##################################################
# fasta indexing                                 #
##################################################
sub index_exists { -e $_.".bwt" }
sub create_index {
  my ($file) = @_;
  my $basename = basename( $file );
  my $index = "$options{'tmp_dir'}/$basename";
  my $cmd = "$options{'bwa_path'} index -p $index $file";
  &run( $cmd );
  return $index;
}
##################################################

##################################################
# Create the bam file                            #
##################################################
sub create_bam {
  my ($ref, $query_files, $prefix) = @_;

  # Remove extension
  my $refname = $1 if( $ref =~ /.*\/([^\/]+)\.[^\.]+$/ );
  die("Could not remove extension from reference filename $ref") unless( $refname );

  my @sais;
  foreach my $in ( @{$query_files} ) {
	my $basename = basename( $in, qw(fastq,fq,txt) );
	my $out = "$options{output_dir}/$refname\_${basename}_aln_sa.sai";
	&aln( $ref, $in, $out, $options_string );
	push(@sais, $out);
  }

  my $basename;
  if ( @sais > 1 ) {
	$basename = &longest_common_prefix( basename($sais[0]), basename($sais[1]) );
  } else {
	$basename = basename( $sais[0], qw(_aln_sa.sai) );
  }
  my $sam_file = "$options{output_dir}/$prefix.sam";
  &sam( $ref, \@sais, $query_files, $sam_file, "-n $options{'num_aligns'}" );
	
  # And convert to bam
  my $basename = basename( $sam_file, '.sam' );
  my $bam_file = "$options{'output_dir'}/$basename.bam";
  &to_bam( $sam_file, $bam_file );

  ## Remove the sam file if the user wants
  &run("rm -f $sam_file") if( $options{'cleanup'} );

  return $bam_file;
}

sub aln {
  my ($ref, $in, $out, $options) = @_;
  my $cmd = "$options{bwa_path} aln $options $ref $in > $out";
  &run($cmd);
}

sub sam {
  my ($ref, $sais, $fastqs, $out, $options) = @_;
  
  ## The number of sais and fastqs should be the same
  die("Expected same number of sai and fastq files. Got ".scalar(@{$sais})." sai files and ".
	  scalar(@{$fastqs})." fastq files") unless( @{$sais} == @{$fastqs} );

  ## Determine if paired end or not
  my $method = (@{$fastqs} > 1) ? "sampe" : "samse";

  my $inputs_string = join(" ", map { "\"$_\"" } (@{$sais}, @{$fastqs}) );
  my $cmd = "$options{bwa_path} $method $options $ref $inputs_string > $out";
  &run($cmd);
}

sub to_bam {
  my ($sam_file, $bam_file) = @_;
  my $cmd = "$options{'samtools_path'} view -bS $sam_file > $bam_file";
  &run($cmd);
}

sub longest_common_prefix {
  my $prefix = shift;
  for( @_ ) {
	chop $prefix while(! /^\Q$prefix\E/);
  }
  $prefix =~ s/_$//;
  return $prefix;
}
##################################################

##################################################
#  Index the bam file                            #
##################################################
# This will replace the existing bam file.
sub index_bam {
  my ($bam) = @_;
  my ($name, $path) = fileparse( $bam, qw(.bam) );
  my $backup = $path."/$name.unsorted.bam";
  &run("mv $bam $backup");
  &bam_sort( $backup, "$path/$name" );
  &bam_index( $bam );
  &run("rm -rf $backup") if( $options{'cleanup'} );
}
sub bam_sort {
  my ($bam, $out) = @_;
  my $cmd = "$options{'samtools_path'} sort $bam $out";
  &run( $cmd );
}

sub bam_index {
  my ($bam) = @_;
  my $cmd = "$options{'samtools_path'} index $bam";
  &run($cmd);
}
##################################################


sub run {
  my ($cmd) = @_;
  system($cmd) == 0 or die("Unable to run $cmd");
}

sub check_parameters {
  
  # display documentation
  if( $options{'help'} ){
    pod2usage( {-exitval=>0, -verbose => 2, -output => \*STDOUT} );
  }
 
  ## make sure we have both tmp and output directories
  foreach my $req ( qw(output_dir tmp_dir ref_file query_file) ) {
	die("Option $req is required") unless( exists( $options{$req} ) );
  }

  # get the query files
  open(IN, "< $options{'query_file'} ") or die("Could not open $options{'query_file'}: $!");
  chomp( my $l = <IN> );
  if ( $l =~ /^\@/ ) {
	push(@{$query_files}, $options{'query_file'});
  } elsif ( -e $l ) {
	push( @{$query_files}, $l );
	chomp( my @therest = <IN> );
	push( @{$query_files}, @therest );
  } else {
	die("Could not determine the format of file: $options{'query_file'}. Expected fastq or list of fastq files");
  }
  close(IN);

  my $query_prefix= $1 if( $options{'query_file'} =~ m|/([^/\.]+)\.[^/]+$| );
  die("Could not parse prefix from $options{'query_file'}") unless( defined( $query_prefix ) );
  
  my $ref_prefix = $1 if( $options{'ref_file'} =~ m|/([^/\.]+)\.[^/]+$| );
  die("Could not parse prefix from $options{'ref_file'}") unless( defined( $ref_prefix ) );

  $prefix = $ref_prefix."_".$query_prefix;

  # only allow 1 or 2 query files
  if ( @{$query_files} < 1 || @{$query_files} > 2 ) {
	die("Found ".scalar(@{$query_files})." query files. Will only accept 1 (single end) or 2 (paired end)");
  }
	
  # The options which will be passed into aln
  my $options_to_param_keys = {
							   'mismatch' => '-M',
							   'max_gaps' => '-o',
							   'max_gap_extensions' => '-e',
							   'open_gap_penalty' => '-O',
							   'extend_gap_penalty' => '-E',
							   'threads' => '-t'
							  };

  my $opts = [];
  foreach my $key (keys %$options_to_param_keys) {
	if ($options{$key}) {
	  push(@$opts, "$options_to_param_keys->{$key} $options{$key}");
	}
  }
  $options_string = join(" ",@$opts);
    
  return 1;
}
