#!/usr/bin/env perl

=head1  NAME 

make-clovr-comparative-circleator-figure.pl - Automatically generate and rasterize a Circleator figure given GenBank format input file(s).

=head1 SYNOPSIS

make-clovr-comparative-circleator-figure.pl
         --circleator_path=/usr/local/packages/circleator/bin/circleator
         --rasterizer_path=/usr/local/packages/circleator/bin/rasterize-svg
         --gb_list_file=genbank.list
         --gene_summary_file=clovr_gene_summary.txt
         --snp_file=somefile_mini.snps
        [--svg_file=figure.svg
         --output_dir=.
         --output_formats='pdf,jpg'
         --output_width=3000
         --output_height=3000
         --help
         --man]

=head1 OPTIONS

B<--gb_list_file,-l>
    Path to the file used in the INPUT_FILE_LIST parameter of the second step (genbank2bsml.default)
    of the CloVR comparative pipeline. The gb_list_file lists the paths to the GenBank files (one
    per line) used as input to the comparative pipeline.

B<--gene_summary_file,-g>
    Path to one of the gene summary files produced by the CloVR comparative pipeline.

B<--snp_file>
    Path to the SNP file produced by the CloVR comparative pipeline.

B<--circleator_path>
    Path to the "circleator" script in the bin/ directory of the Circleator distribution.

B<--rasterizer_path>
    Path to the "rasterize-svg" script in the bin/ directory of the Circleator distribution.

B<--svg_file,-s>
    Optional. Name of the SVG file to generate in --output_dir. Default is "figure.svg"

B<--output_dir>
    Optional. Path to the directory in which the contig list file, Circleator configuration file,
    output SVG file, and rasterized file(s) should be written. Default is current working directory.

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

This is a wrapper script for Circleator that performs the following tasks:

1. Parses the CloVR-comparative gene summary file specified by --gene_summary_file.
2. Identifies the set of GenBank files that correspond to the reference genes in --gene_summary_file
3. Creates a Circleator --contig_list file for all of the GenBank files identified in the previous step.
4. Converts the gene cluster information in --gene_summary_file into a standard Circleator format.
5. Creates a custom Circleator configuration file that displays:
  a. The reference genes (forward and reverse strand)
  b. The "core" genes (those present in all strains according to the gene clusters)
  c. The unique genes (those present only in the reference strain according to the gene clusters)
  d. The unique SNPs (positions where the reference differs from _all_ of the other strains)
6. Runs the Circleator to generate an SVG-format figure
7. Runs the SVG rasterizer to convert the SVG file into PDF, PNG, or JPEG

=head1 INPUT

The following input and output files from an instance of the CloVR-Comparative ergatis pipeline:
 o the INPUT_FILE_LIST of the genbank2bsml.default component
   i.e., a list of GenBank files
 o one of the gene summary files produced by the pipeline
 o the SNP output file produced by the pipeline

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

# If the SNP file contains more than this number of SNPs then the SNPs will not
# be included in the resulting figure. This is a workaround that is unlikely
# to be triggered now that the SNP files are being pre-filtered for unique SNPs.
my $MAX_SNP_FILE_SIZE = 150000;

## input
my $options = {};

&GetOptions($options, 
            "gb_list_file|l=s",
            "gene_summary_file|g=s",
            "snp_file=s",
            "circleator_path=s",
            "rasterizer_path=s",
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
my($gb_list_file, $gene_summary_file, $snp_file) = map {$options->{$_}} ('gb_list_file', 'gene_summary_file', 'snp_file');

# read gb_list_file, parse contig ids from genbank files
my %contig2gbacc = (); # adding hash of genbank accession id to gbk files because sometimes accession id is used in gene summary files instead of locus
my $contig2gb = &read_gb_list_file($gb_list_file);

my $fh = FileHandle->new();
$fh->open($gene_summary_file) || die "unable to read from CloVR-Comparative gene summary file $gene_summary_file";
my $lnum = 0;
my $num_data_cols = undef;
my $num_strains = undef;
my $strain_names = [];
my $ref_strain = undef;

# locations of cluster and SNP data, respectively
my $first_cluster_col = undef;
my $last_cluster_col = undef;
my $first_snp_col = undef;
my $last_snp_col = undef;
my $expected_cols = undef;

# track GenBank files in which reference genes appear
my $ref_gb_files = {};

# gene clusters in the order they appear in the file
my $clusters = [];
# gene clusters indexed by ref gene - used to find duplicates
my $rg2cluster = {};

while (my $line = <$fh>) {
  my @f = split(/\t/, $line);
  chomp($f[-1]);
  my $nf = scalar(@f);
  ++$lnum;
  chomp($line);
  
  if ($lnum == 1) {
    if ($line =~ /^Gene Info\t{6}Genes In Cluster(\t+)SNPs\s*$/) {
      $num_strains = length($1);
      $first_cluster_col = 6;
      $last_cluster_col = 6 + $num_strains - 1;
      $first_snp_col = $last_cluster_col + 1;
      $last_snp_col = $first_snp_col + $num_strains - 1;
      $expected_cols = $last_snp_col + 1;
    } else {
      die "unexpected value at line $lnum: $line";
    }
  } 
  # parse strain/genome names from secondary column headers
  elsif ($lnum == 2) {
    $num_data_cols = $nf;
    die "unexpected number of columns ($nf instead of $expected_cols) at line $lnum" unless ($nf == $expected_cols);
    for (my $c = $first_cluster_col;$c <= $last_cluster_col;++$c) {
      push(@$strain_names, $f[$c]);
    }
    my $ncg = 0;
    # use "Current Genome" SNP header to determine the current genome/strain
    my $strains_h = {};
    map {$strains_h->{$_} = 1;} @$strain_names;
    for (my $c = $first_snp_col;$c <= $last_snp_col;++$c) {
      if ($f[$c] ne 'Current Genome') {
        delete $strains_h->{$f[$c]};
      } else {
        ++$ncg;
      }
    }
    die "found $ncg 'Current Genome' column(s), expected 1" if ($ncg != 1);

    my @rk = keys(%$strains_h);
    if (scalar(@rk) == 1) {
      $ref_strain = $rk[0];
    } else {
      die "unable to determine reference strain from SNP headers";
    }

    die "no reference strain defined" if (!defined($ref_strain));
  } 
  # data
  else {
    my $has_snps = 1;
    my $has_cluster = 1;
    my $strain2genes = {};

    # no SNP columns
    if ($nf == ($num_data_cols - $num_strains)) {
      $has_snps = 0;
    } 
    # unexpected column count
    elsif ($nf != $num_data_cols) {
      die "unexpected number of columns ($nf instead of $num_data_cols) at line $lnum" unless ($nf == $num_data_cols);
    }
    # lookup GenBank file for reference gene/contig
    my $ref_contig = $f[0];
    my $ref_gb = $contig2gb->{$ref_contig};
    if (!defined($ref_gb)) {
	if(defined($contig2gbacc{$ref_contig})) {
		$ref_gb = $contig2gbacc{$ref_contig};
	} else {
    		print STDERR "ERROR - reference contig $ref_contig not found in any GenBank file listed in $gb_list_file\n";
	}
    }
    ++$ref_gb_files->{$ref_gb};

    # parse cluster info
    my $ng = 0;
    for (my $c = $first_cluster_col;$c <= $last_cluster_col;++$c) {
      my $cv = $f[$c];
      my $sn = $strain_names->[$c - $first_cluster_col];
      my $list = [];
      # at least one gene
      if ($cv !~ /^-$/) {
        # delimiters: 
        #  -comma separates lists of genes grouped by molecule
        #  -colon separates molecule name from list of genes
        #  -genes are separated by "|"
        #  -each gene consists of molecule_name ||| gene_name (an unfortunate choice given the earlier use of "|")

        # replace within-gene delimiters
        $cv =~ s/\|\|\|/!!!/g;
        my @mols = split(/\,/, $cv);
        foreach my $mol (@mols) {
          my($moln, $gene_list) = ($mol =~ /^([^:]+):(.*)$/);
          my @genes = split(/\|/, $gene_list);
          foreach my $gene (@genes) {
            my($mn, $gn) = split(/\!{3}/, $gene);
            push(@$list, $gn);
          }
        }
      }
      $strain2genes->{$sn} = $list;
    }
    # include ref gene to correct strain
    my $ref_list = $strain2genes->{$ref_strain};
    die "couldn't find reference strain '$ref_strain'" if (!defined($ref_list));
    my($rm, $rg) = split(/\|\|\|/, $f[1]);
    push(@$ref_list, $rg);

    foreach my $strain (keys %$strain2genes) {
      my $list = $strain2genes->{$strain};
      my $nsg = scalar(@$list);
      $ng += $nsg;
      my $isref = ($strain eq $ref_strain) ? " (REF)" : "";
    }

    # merge clusters with the same reference gene
    #  (since it looks like the clusters are currently computed per-fragment, not per-gene)
    my $c = $rg2cluster->{$rg};
    if (defined($c)) {
      print STDERR "INFO - merging cluster with identical reference gene $rg at line $lnum\n";
      push(@{$c->{'ref_molecules'}}, $rm);
      foreach my $key (keys %$strain2genes) {
        my $sl = $strain2genes->{$key};
        my $ol = $c->{'strain2genes'}->{$key};
        # remove duplicates
        my $nh = {};
        map { $nh->{$_} = 1; } (@$sl, @$ol);
        my @new_list = keys %$nh;
        $c->{'strain2genes'}->{$key} = \@new_list;
      }
    } 
    # new cluster
    else {
       $c = { 'ref_molecules' => [$rm], 'ref_gene' => $rg, 'strain2genes' =>  $strain2genes };
      push(@$clusters, $c);
      $rg2cluster->{$rg} = $c;
    }
  }
}
$fh->close();

# $ref_gb_files contains the list of reference files to pass to Circleator
my @ref_gb_files = keys %$ref_gb_files;
foreach my $gbf (@ref_gb_files) {
  my $count = $ref_gb_files->{$gbf};
  print "INFO - found $count ref gene(s) in " . $gbf . "\n";
}

# create Circleator contig list file using same methodology as CloVR-Microbe version
# read length of each sequence in each .gbf file (assumes 1 sequence per file)
my $file2seqlen = {};
my $total_seqlen = 0;
# also may need to generate unique seq ids
my $file2seqid = {};
my $seqids = {};
foreach my $gbf_file (@ref_gb_files) {
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
my $conf_file = $svg_file;
$conf_file =~ s/\.svg$/.conf/;
my $cluster_file = $svg_file;
$cluster_file =~ s/\.svg$/-gene-clusters.txt/;

my $svg_file_path = File::Spec->catfile($output_dir, $svg_file);
my $list_file_path = File::Spec->catfile($output_dir, $list_file);
my $log_file_path = File::Spec->catfile($output_dir, $log_file);
my $conf_file_path = File::Spec->catfile($output_dir, $conf_file);
my $cluster_file_path = File::Spec->catfile($output_dir, $cluster_file);

# generate contig list file 
&write_contig_list_file($list_file_path, \@sorted_files, $file2seqlen);

# create Circleator config file
&write_config_file($conf_file_path, $cluster_file_path, $strain_names, $ref_strain, $snp_file);

# write Cluster_Table fileconf
&write_cluster_file($cluster_file_path, $strain_names, $ref_strain, $clusters);

# run Circleator (param)
my $circleator = $options->{'circleator_path'};
my $c_cmd = "$circleator --pad=550 --contig_gap_size=5000 --config=$conf_file_path --contig_list=$list_file_path --log=$log_file_path --debug=all >$svg_file_path";
&run_sys_command($c_cmd);

# run rasterizer once for each format in --output_formats
my $rasterizer = $options->{'rasterizer_path'};
my($ofs, $ow, $oh) = map {$options->{$_}} ('output_formats', 'output_width', 'output_height');
my @output_formats = split(/\s*\,\s*/, $ofs);

foreach my $of (map {lc($_)} @output_formats) {
  if ($of !~ /^pdf|png|jpg$/) {
    print STDERR "WARN - skipping unrecognized output format $of\n";
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
  my ($file_base,$file_dir,$file_ext);
  ## make sure required parameters were passed
  my @required = qw(gb_list_file gene_summary_file snp_file circleator_path rasterizer_path);
  for my $option ( @required ) {
    unless ( defined $options->{$option} ) {
      die("--$option is a required option");
    }
  }

  ## defaults
  $options->{'output_formats'}= $DEFAULT_OUTPUT_FORMATS if (!defined($options->{'output_formats'}));
  $options->{'output_width'}= $DEFAULT_OUTPUT_WIDTH if (!defined($options->{'output_width'}));
  $options->{'output_height'}= $DEFAULT_OUTPUT_HEIGHT if (!defined($options->{'output_height'}));
  ## Added code to assign organism name as filename instead of prefix from clovr comparative pipeline
  if (!defined($options->{'svg_file'})) {
  	$options->{'svg_file'}= $DEFAULT_SVG_FILE;
  } else {
  	($file_base,$file_dir,$file_ext) = fileparse($options->{'svg_file'},qr/\.[^.]*/);
	$options->{'svg_file'} = $file_base.".svg";	
  }
  $options->{'output_dir'}= $DEFAULT_OUTPUT_DIR if (!defined($options->{'output_dir'}));

  ## additional parameter checking
}

sub get_file_line_count {
  my($file) = @_;
  my $count = 0;
  my $fh = FileHandle->new();
  $fh->open($file) || die "unable to read from $file";
  while (my $line = <$fh>) {
    ++$count;
  }
  $fh->close();
  return $count;
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
  my $nf = scalar(@$sorted_files);
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

sub read_gb_list_file {
  my($file) = @_;
  my $gb_files = [];

  my $fh = FileHandle->new();
  $fh->open($file) || die "unable to read from gb_list file $file";
  my $lnum = 0;
  my $nnf = 0;
  while (my $line = <$fh>) {
    chomp($line);
    if (!(-e $line) || !(-r $line)) {
      print STDERR "ERROR - $line does not exist or is not readable\n";
      ++$nnf;
    } else {
      push(@$gb_files, $line);
    }
  }
  $fh->close();
  die "one or more GenBank files not found" if ($nnf > 0);

  # contig 2 GenBank file mapping
  my $contig2gb = {};
  foreach my $gbf (@$gb_files) {
    my $contig_ids_str = `egrep '^LOCUS' $gbf`;
    my $contig_alt_id = `egrep '^ACCESSION' $gbf`;
    foreach my $locus_line (split(/\n/, $contig_ids_str)) {
      if ($locus_line =~ /^LOCUS\s+(\S+)\s+\d+ bp/) {
        $contig2gb->{$1} = $gbf;
      } 
      # error - contig id is too long and collides with sequence length
      elsif ($locus_line =~ /^LOCUS\s+(\S+)\d+ bp/) {
        print STDERR "ERROR - malformed LOCUS line in $file: $locus_line\n";
      }
      else {
        print STDERR "ERROR - unable to parse LOCUS line in $file: $locus_line\n";
      }
    }
    foreach my $acc_line (split(/\n/, $contig_alt_id)) {
    	if($acc_line =~ /^ACCESSION\s+(\S+)/) {
		$contig2gbacc{$1} = $gbf;
    	}
    }
  }
  return $contig2gb;
}

sub write_cluster_file {
  my($cluster_file_path, $strain_names, $ref_strain, $clusters) = @_;
  my $fh = FileHandle->new();
  $fh->open(">$cluster_file_path") || die "unable to write to gene cluster file at $cluster_file_path";

  # print header line
  $fh->print(join("\t", 'ref_gene', @$strain_names) . "\n");

  # print gene clusters
  foreach my $cluster (@$clusters) {
    my($rg, $s2g) = map {$cluster->{$_}} ('ref_gene', 'strain2genes');
    $fh->print($rg . "\t");
    my $sl = [];
    foreach my $sn (@$strain_names) {
      my $list = $s2g->{$sn};
      $list = [] if (!defined($list));
      push(@$sl, join(",", @$list));
    }
    $fh->print(join("\t", @$sl));
    $fh->print("\n");
  }
  $fh->close();
}

sub write_config_file {
  my($config_file_path, $cluster_file, $strain_names, $ref_strain, $snp_file) = @_;
  my $strain_str = join("|", @$strain_names);
  my $core_sig = "";
  my $unique_sig = "";

  foreach my $strain (@$strain_names) {
    $core_sig .= "1";
    $unique_sig .= ($strain eq $ref_strain) ? "1" : "0";
  }

  # num strains minus 1
  my $nsm1 = scalar(@$strain_names) - 1;

  # only display SNPs if SNP line count is < $MAX_SNP_FILE_SIZE
  my $num_snps = &get_file_line_count($snp_file) - 1;
  my $snp_block = "";
  if ($num_snps <= $MAX_SNP_FILE_SIZE ) {
    $snp_block = <<SNP_BLOCK;
# load SNPs
new ls1 load feat-file=${snp_file},feat-file-type=snp-table,snp-ref=$ref_strain
# display SNPs unique to the reference
# only display SNPs where all the other genomes differ from the reference:
new us1 rectangle heightf=0.06,feat-type=SNP,color1=#4030ff,color2=#4030ff,feat-tag=SNP_num_diffs,feat-tag-value=$nsm1
small-label label-text=unique&nbsp;SNPs
SNP_BLOCK
  }

  my $fh = FileHandle->new();
  $fh->open(">$config_file_path") || die "unable to write to configuration file at $config_file_path";
  print $fh <<CONFIG;
# Circleator config file for CloVR-Comparative

new lgc1 load-gene-cluster-table gene-cluster-file=$cluster_file

coords
tiny-cgap

contigs c1
small-label label-text=contigs
tiny-cgap

new cds-fwd rectangle 0.06 . . . CDS 1 #000000
rRNAs-fwd color2=#00ff00,stroke-width=1,innerf=same,outerf=same
tRNAs-fwd color2=#ff0000,stroke-width=1,innerf=same,outerf=same
tiny-cgap
new cds-rev rectangle 0.06 . . . CDS -1 #000000
rRNAs-rev color2=#00ff00,stroke-width=1,innerf=same,outerf=same
tRNAs-rev color2=#ff0000,stroke-width=1,innerf=same,outerf=same
small-label label-text=genes&nbsp;and&nbsp;rRNAs(green)&nbsp;and&nbsp;tRNAs(red)
small-cgap

# core genes (conserved in all compared genomes)
new cg1 rectangle 0.06 feat-type=gene,gene-cluster-genomes=$strain_str,gene-cluster-signature=$core_sig,color1=#d230ff,color2=#d230ff
small-label label-text=core&nbsp;genes
small-cgap

# unique genes (conserved only in reference genome)
new ug1 rectangle 0.06 feat-type=gene,gene-cluster-genomes=$strain_str,gene-cluster-signature=$unique_sig,color1=#d2612a,color2=#d2612a
small-label label-text=unique&nbsp;genes
small-cgap

$snp_block

medium-cgap

%GCmin-max color2=#ff0000
small-label label-text=%GC
tiny-cgap
GCskew-min-max-df0 color2=#00ff00
small-label label-text=GC-skew

# label only contigs that are at least 50kb
medium-label heightf=0.04,feat-track=c1,label-function=primary_id,packer=none,label-type=spoke,innerf=1.1,feat-min-length=20000,feat-type=contig

CONFIG
  $fh->close();
}

# DEBUG - expand a single contig
#new ssl1 scaled-segment-list user-feat-fmin=0,user-feat-fmax=50000,user-feat-type=roi,scale=140

# DEBUG
#new cg1 rectangle 0.06 feat-type=gene,gene-cluster-genomes=Escherichia_coli_3006|Escherichia_coli_34870|Escherichia_coli_52239|Escherichia_coli_5412,gene-cluster-signature=1001,color1=#0000ff

# DEBUG - highlight and label genes in expanded region
#new egl1 label overlapping-feat-type=roi,feat-type=gene,packer=none,label-function=locus,heightf=0.03,innerf=1,label-type=spoke
#new eglh1 rectangle overlapping-feat-type=roi,feat-type=gene,innerf=0.1,outerf=1.25,opacity=0.1,color1=#00ff00,color2=#000000
