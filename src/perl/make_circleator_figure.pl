#!/usr/bin/env perl

=head1  NAME 

make-circleator-figure.pl - Automatically generate and rasterize a Circleator figure given GenBank format input file(s).

=head1 SYNOPSIS

make-circleator-figure.pl
         --genbank_dir=.
         --circleator_path=/usr/local/packages/circleator/bin/circleator
         --rasterizer_path=/usr/local/packages/circleator/bin/rasterize-svg
         --config_dir=/path/to/predefined/config/file/directory
        [--svg_file=figure.svg
         --output_dir=.
         --output_formats='pdf,jpg'
         --output_width=3000
         --output_height=3000
         --help
         --man]

=head1 OPTIONS
    
B<--genbank_dir,-g>
    Path to a directory that contains one or more GenBank flat files with the ".gbf" suffix.
    The named directory will be searched recursively for all such files.

B<--circleator_path>
    Path to the "circleator" script in the bin/ directory of the Circleator distribution.

B<--rasterizer_path>
    Path to the "rasterize-svg" script in the bin/ directory of the Circleator distribution.

B<--config_dir,-c>
    Path to a directory that contains the following predefined Circleator config files:
     -fig-1.cfg
     -fig-1-small.cfg

B<--svg_file,-s>
    Optional. Name of the SVG file to generate in --output_dir. Default is "figure.svg"

B<--output_dir>
    Optional. Path to the directory in which the contig list file, output SVG file, and
    rasterized file should be written. Default is current working directory.

B<--output_formats>
    Optional. Comma-delimited list of formats ('png', 'jpg', or 'pdf') to which SVG should 
    be rasterized/converted. Default is 'pdf,jpg'

B<--output_width>
    Optional. Width of the final rasterized/converted images. Default is 3000.

B<--output_height>
    Optional. Height of the final rasterized/converted images. Default is 3000.

B<--help,-h> 
    Display the script documentation.

B<--man,-m>
    Display the script documentation.

=head1 DESCRIPTION

Wrapper script for Circleator that performs the following tasks:

1. Creates a Circleator --contig_list file for all the .gbf output files in a given directory.
2. Selects one of 3 predefined config files based on the number of input seqs and their lengths
3. Runs the Circleator to generate an SVG-format figure
4. Runs the SVG rasterizer to convert the SVG file into PDF, PNG, or JPEG

=head1 INPUT

One or more GenBank-format annotation files.

=head1 OUTPUT

A graphical plot of the input sequence(s) in both SVG format and whatever format is specified
by the --output_format option.

=head1 CONTACT

    Jonathan Crabtree
    jonathancrabtree@gmail.com

=cut

use strict;
use File::Basename;
use FileHandle;
use File::Spec;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Pod::Usage;

## globals
my $DEFAULT_OUTPUT_FORMATS = 'pdf,jpg';
my $DEFAULT_OUTPUT_WIDTH = 3000;
my $DEFAULT_OUTPUT_HEIGHT = 3000;
my $DEFAULT_SVG_FILE = 'figure.svg';
my $DEFAULT_OUTPUT_DIR = '.';

## input
my $options = {};

&GetOptions($options, 
            "genbank_dir|g=s",
            "circleator_path=s",
            "rasterizer_path=s",
            "config_dir|c=s",
            "svg_file|s=s",
            "output_dir=s",
            "output_formats=s",
            "output_width=i",
            "output_height=i",
            "help|h",
            "man|m"
           );

## display documentation
if ( $options->{'help'} || $options->{'man'} ) {
  pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDERR} );
}
&check_parameters($options);

## main program

# traverse genbank_dir to find .gbf files
my $genbank_dir = $options->{'genbank_dir'};
my $gbf_files_str = `find $genbank_dir -type f -name '*.gbf' -print`;
my @gbf_files = split(/\n/, $gbf_files_str);
my $nf = scalar(@gbf_files);
print STDERR "found $nf .gbf file(s) in $genbank_dir\n";

# read length of each sequence in each .gbf file (assumes 1 sequence per file)
my $file2seqlen = {};
my $total_seqlen = 0;
# also may need to generate unique seq ids
my $file2seqid = {};
my $seqids = {};
foreach my $gbf_file (@gbf_files) {
  my ($seqlen, $seqid) = &get_genbank_seqlen_plus_accession($gbf_file);
#  print STDERR "$gbf_file seqlen=$seqlen seqid=$seqid\n";
  $file2seqlen->{$gbf_file} = $seqlen;
  $total_seqlen += $seqlen;
  if (defined($seqids->{$seqid})) {
    ++$seqids->{$seqid};
    $seqid = "unknown" . sprintf("%05d", $seqids->{$seqid});
  }
  $file2seqid->{$gbf_file} = $seqid;
  $seqids->{$seqid} = 1;
}

# sort by decreasing length
my @sorted_files = sort { $file2seqlen->{$b} <=> $file2seqlen->{$a}} keys %$file2seqlen;

# output file paths
my $output_dir = $options->{'output_dir'};
my $svg_file = $options->{'svg_file'};
my $list_file = $svg_file;
$list_file =~ s/\.svg$/-contig-list.txt/;
my $log_file = $svg_file;
$log_file =~ s/\.svg$/.log/;
my $svg_file_path = File::Spec->catfile($output_dir, $svg_file);
my $list_file_path = File::Spec->catfile($output_dir, $list_file);
my $log_file_path = File::Spec->catfile($output_dir, $log_file);

# generate contig list file 
&write_contig_list_file($list_file_path, \@sorted_files, $file2seqlen);

# select config file based on number of contigs and total sequence length
my $conf_file = undef;

# single scaffold/contig
if ($nf == 1) {
  if ($total_seqlen <= 500000) {
    $conf_file = "fig-1-small.cfg";
  } else {
    $conf_file = "fig-1.cfg";
  }
}
# multiple scaffolds/contigs
else {
  $conf_file = "fig-1-multiple-contigs.cfg";
}

my $conf_file_path = File::Spec->catfile($options->{'config_dir'}, $conf_file);

# run Circleator (param)
my $circleator = $options->{'circleator_path'};
my $c_cmd = "$circleator --pad=500 --config=$conf_file_path --contig_list=$list_file_path --log=$log_file_path >$svg_file_path";
&run_sys_command($c_cmd);

# run rasterizer once for each format in --output_formats
my $rasterizer = $options->{'rasterizer_path'};
my($ofs, $ow, $oh) = map {$options->{$_}} ('output_formats', 'output_width', 'output_height');
my @output_formats = split(/\s*\,\s*/, $ofs);

foreach my $of (map {lc($_)} @output_formats) {
  if ($of !~ /^pdf|png|jpg$/) {
    print STDERR "skipping unrecognized output format $of\n";
    next;
  }
  # special case for PDF format: using large width/height values can introduce artifacts
  my($w, $h) = ($of eq 'pdf') ? (1000,1000) : ($ow, $oh);
  my $r_cmd = "$rasterizer $svg_file_path $of $w $h";
  &run_sys_command($r_cmd);
}

exit(0);

## subroutines

sub check_parameters {
  my $options = shift;
    
  ## make sure required parameters were passed
  my @required = qw(genbank_dir circleator_path rasterizer_path config_dir);
  for my $option ( @required ) {
    unless ( defined $options->{$option} ) {
      die("--$option is a required option");
    }
  }

  ## defaults
  $options->{'output_formats'}= $DEFAULT_OUTPUT_FORMATS if (!defined($options->{'output_formats'}));
  $options->{'output_width'}= $DEFAULT_OUTPUT_WIDTH if (!defined($options->{'output_width'}));
  $options->{'output_height'}= $DEFAULT_OUTPUT_HEIGHT if (!defined($options->{'output_height'}));
  $options->{'svg_file'}= $DEFAULT_SVG_FILE if (!defined($options->{'svg_file'}));
  $options->{'output_dir'}= $DEFAULT_OUTPUT_DIR if (!defined($options->{'output_dir'}));

  ## additional parameter checking
}

# Get length and name of the first sequence in a GenBank flat file.
sub get_genbank_seqlen_plus_accession {
  my($file) = @_;
  my $seqlen = undef;
  my $accession = undef;
  my $fh = FileHandle->new();
  $fh->open($file) || die "unable to read from GenBank flat file $file";
  while (my $line = <$fh>) {
    if ($line =~ /^LOCUS\s+(\S+)\s+(\d+) bp/) {
      $accession = $1;
      $seqlen = $2;
    }
    elsif (!defined($accession) && ($line =~ /^ACCESSION\s+(\S+)/)) {
      $accession = $1;
    }
    # first line may be unreliable due to long seq ids overrunning seqlen, so look for 'source' specification
    elsif ($line =~ /^     source          1..(\d+)/) {
      $seqlen = $1 if (!defined($seqlen));
      last;
    }
  }
  $fh->close();
  # use filename as accession if one can't be parsed from the file
  if (!defined($accession)) {
    my $file_base = basename($file);
    # special case for CloVR Microbe pipeline
    if ($file_base =~ /^(NODE\_\d+)\_length\_\d+\_cov\_\d+\.(\d+)\.gbf/) {
      $accession = $1. '.' . $2;
    } elsif (/^(.+)\.gbf/) {
      $accession = $1;
    }
    
  }
  return ($seqlen, $accession);
}

sub write_contig_list_file {
  my ($list_file_path, $sorted_files, $file2seqlen) = @_;
  my $cfh = FileHandle->new();
  $cfh->open(">$list_file_path") || die "unable to write to $list_file_path";

  my $get_gap_size = sub {
    my($seqlen) = @_;
    # 5kb or 0.25 * the contig length, whichever is smaller
    my $hcl = int($seqlen * 0.25);
    return ($hcl < 5000) ? $hcl : 5000;
  };

  for (my $i = 0;$i < $nf;++$i) {
    my $sf = $sorted_files->[$i];
    my $seqid = $file2seqid->{$sf};
    my $seqlen = $file2seqlen->{$sf};
    $cfh->print(join("\t", ($seqid, "", "", $sf, "", "")). "\n");
    # insert gap after each contig in multi-contig figure
    if ($nf > 1) {
      # gap size is based on the length of the current contig and the next one
      my $next = ($i < ($nf - 1)) ? $i + 1 : 0;
      my $next_seqlen = $file2seqlen->{$sorted_files->[$next]};
      my $gap_size = &$get_gap_size($seqlen) + &$get_gap_size($next_seqlen);
      $cfh->print(join("\t", ("gap", "", $gap_size, "", "", "")). "\n");
    }
  }
  $cfh->close();
}

sub run_sys_command {
  my($cmd) = @_;
  system($cmd);

  # check for errors, halt if any are found
  my $err = undef;
  if ($? == -1) {
    $err = "failed to execute: $!";
  }
  elsif ($? & 127) {
    $err = sprintf("child died with signal %d, %s coredump\n", ($? & 127), ($? & 128) ? 'with' : 'without');
  }
  else {
    my $exit_val = $? >> 8;
    $err = sprintf("child exited with value %d\n", $exit_val) if ($exit_val != 0);
  }
  die $err if (defined($err));
}

