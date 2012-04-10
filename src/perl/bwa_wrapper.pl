#!/usr/bin/perl

=head1  NAME 

lgt_bwa.pl

=head1 SYNOPSIS


      
=head1 OPTIONS

=over 8

This help message

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
						  'query_list|ql=s',
						  'ref_file|r=s',
						  'ref_file_list|rl=s',
						  'output_dir|o=s',
						  'mismatch|mm=s',
						  'max_gaps|mg=s',
						  'max_gap_extensions|mge=s',
						  'open_gap_penalty|og=s',
						  'extend_gap_penalty|eg=s',
						  'threads|t=s',
						  'use_bwasw|ub=s',
						  'num_aligns|na=s',
						  'bwa_path|b=s',
						  'samtools_path|st=s',
						  'cleanup|c=s',
						  'help|h');

# display documentation
if( $options{'help'} ){
    pod2usage( {-exitval=>0, -verbose => 2, -output => \*STDOUT} );
}

my $PAIRED = 0;
my $ref_files = [];
my $query_files = [];
my $options_string = '';

## make sure all passed options are peachy
&check_parameters(\%options);

if($options{use_bwasw}) {
    die "bwasw is not implemented yet\n";
}
else {
    &run_bwa();
}

sub run_bwa {

  foreach my $ref (@$ref_files) {
	#Is the reference indexed? If not, index it.
	my $temp_index = 0;
	unless( index_exists( $ref ) ) {
	  $ref = &create_index( $ref );
	  $temp_index = 1;
	}

	chomp $ref;
	$ref =~ /.*\/([^\/]+)\.[^\.]+$/;
	my $refname = $1;

	my $sam_file;

	# In here if we're paired
	if ($PAIRED) {

	  my $in1 = $query_files->[0];
	  my $basename1 = basename( $in1, qw(fastq,fq,txt) );
	  my $out1 = "$options{output_dir}/${refname}_${basename1}_aln_sa1.sai";
	  my $in2 = $query_files->[1];
	  my $basename2 = basename( $in2, qw(fastq,fq,txt) );
	  my $out2 = "$options{output_dir}/${refname}_${basename2}_aln_sa2.sai";

            
	  # Run the first one through aln
	  my $cmd = "$options{bwa_path} aln $options_string $ref $in1 > $out1";
	  print "Running: $cmd\n";
	  system($cmd) == 0 or die "Unable to run $cmd\n";
            
	  # Run the second one through aln
	  $cmd = "$options{bwa_path} aln $options_string $ref $in2 > $out2";
	  print "Running: $cmd\n";
	  system($cmd) == 0 or die "Unable to run $cmd\n";

	  # Run sampe
	  my $prefix = &longest_common_prefix( basename($in1), basename($in2) );
	  $prefix =~ s/_$//;
	  $sam_file = "$options{output_dir}/$refname\_$prefix.sam";
	  $cmd = "$options{bwa_path} sampe -n $options{num_aligns} $ref \"$out1\" \"$out2\" \"$in1\" \"$in2\" > $sam_file";
	  print "Running: $cmd\n";
	  system($cmd) == 0 or die "Unable to run $cmd\n";

	  if ($options{cleanup}) {
		$cmd = "rm -f $out1 $out2";
		print "Running: $cmd\n";
		system($cmd) == 0 or die "Unable to run $cmd\n";
	  }

	} else {

	  my $in = $query_files->[0];
	  my $basename = basename( $in, qw(fastq,fq,txt) );
	  my $out = "$options{output_dir}/$refname\_${basename}_aln_sa.sai";
	  my $cmd = "$options{bwa_path} aln $options_string $ref $in > $out";

	  print "Running: $cmd\n";
	  system($cmd) == 0 or die "Unable to run $cmd\n";

	  $sam_file = "$options{output_dir}/$refname\_${basename}.sam";
	  $cmd = "$options{bwa_path} samse -n $options{num_aligns} $ref \"$out\" \"$in\" > $sam_file";
	  print "Running: $cmd\n";
	  system($cmd) == 0 or die "Unable to run $cmd\n";

	  if ($options{cleanup}) {
		$cmd = "rm -f $out";
		print "Running: $cmd\n";
		system($cmd) == 0 or die "Unable to run $cmd\n";
	  }
	}
	
	# And convert to bam
	my $basename = basename( $sam_file, '.sam' );
	my $bam_file = "$options{'output_dir'}/$basename.bam";
	my $cmd = "$options{'samtools_path'} view -bS $sam_file > $bam_file";
	print "Running: $cmd\n";
	system($cmd) == 0 or die("Unable to run $cmd");

	if( $options{'cleanup'} ) {
	  $cmd = "rm -f $sam_file";
	  print "Running: $cmd\n";
	  system($cmd) == 0 or die "Unable to run $cmd\n";

	  if( $temp_index ) {
		$cmd = "rm -f $ref*";
		print "Running: $cmd\n";
		system($cmd) == 0 or die "Unable to run $cmd\n";
	  }
	}
  }

}

sub index_exists { -e $_.".bwt" }
sub create_index {
  my ($file) = @_;
  my $basename = basename( $file );
  my $index = "$options{'output_dir'}/$basename";
  my $cmd = "$options{'bwa_path'} index -p $index $file";
  system($cmd) && die("Could not run cmd $cmd: $!");
  return $index;
}

sub longest_common_prefix {
  my $prefix = shift;
  for( @_ ) {
	chop $prefix while(! /^\Q$prefix\E/);
  }
  $prefix;
}

sub check_parameters {
    
    if($options{ref_file}) {
        my @files = split(/,/,$options{ref_file});
        $ref_files = \@files;
    }
    elsif($options{ref_file_list}) {
        my @lines = `cat $options{ref_file_list}`;
        $ref_files = \@lines;
    }
    else {
        die "No reference file specified in ref_file or ref_file_list"
    }

	my @qfiles;
	if( $options{'query_file'} ) {
	  @qfiles = split(/,/,$options{'query_file'});
	} else {
	  die("No query file specified.");
	}

	foreach my $qf ( @qfiles ) {
	  open(IN, "< $qf") or die("Could not open $qf: $!");
	  chomp( my $l = <IN> );
	  if( $l =~ /^\@/ ) {
		push(@{$query_files}, $qf);
	  } elsif( -e $l ) {
		push( @qfiles, $l );
		chomp( my @therest = <IN> );
		push( @qfiles, @therest );
	  } else {
		die("Could not determine the format of file: $qf");
	  }
	  close(IN);
	}

	# only allow 1 or 2 query files
	if( @{$query_files} == 2 ) {
	  $PAIRED=1;
	} elsif( @{$query_files} != 1 ) {
	  die("Found ".scalar(@{$query_files})." query files. Will only accept 1 (single end) or 2 (paired end)");
	}
	  

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
        if($options{$key}) {
            push(@$opts, "$options_to_param_keys->{$key} $options{$key}");
        }
    }
    $options_string = join(" ",@$opts);
    
    return 1;
}
