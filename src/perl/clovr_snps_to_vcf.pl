#!/usr/bin/perl

=head1  NAME 

clovr-snps-to-vcf.pl - Convert CloVR-Comparative SNP output file to a set of VCF 4.2.

=head1 SYNOPSIS

clovr-snps-to-vcf.pl
         --snp_file=somefile_mini.snps
        [--output_dir=./
         --help
         --man]

=head1 OPTIONS

B<--snp_file>
    Path to the SNP file produced by the CloVR comparative pipeline.

B<--output_dir>
    Directory to which VCF files should be written. Default is current working directory.

B<--help,-h> 
    Display the script documentation.

B<--man,-m>
    Display the script documentation.

=head1 DESCRIPTION

Script to convert a CloVR-Comparative SNP output file to VCF 4.2.

=head1 INPUT

The SNP output file produced by an instance of the CloVR-Comparative ergatis pipeline.

=head1 OUTPUT

A corresponding VCF file.

=head1 CONTACT

    Jonathan Crabtree
    jonathancrabtree@gmail.com

=cut

use strict;
use FileHandle;
use File::Spec;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Pod::Usage;

## globals
my $VCF_COLS = [ "#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT" ];

## input
my $options = {};

&GetOptions($options, 
            "snp_file=s",
            "output_dir=s",
            "help|h",
            "man|m"
           );

## display documentation
if ( $options->{'help'} || $options->{'man'} ) {
  pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDERR} );
}
&check_parameters($options);

## main program
my($snp_file, $odir) = map {$options->{$_}} ('snp_file', 'output_dir');

# parse the CloVR-Comparative SNP file
my $snp_data = &read_clovr_snps($snp_file);

# loop over reference sequences, writing a VCF file for each, until there are no SNPs left
my($ns, $strains, $snps) = map { $snp_data->{$_} } ('n_strains', 'strains', 'snps');
my $n_indels = 0;

# count number of SNPs mapped to each reference
my @snp_counts = map { 0 } 1 .. $ns;
foreach my $snp (@$snps) {
  my($is_indel, $bases) = map {$snp->{$_}} ('is_indel', 'bases');
  if ($is_indel) {
    ++$n_indels;
    next;
  }

  for (my $s = 0;$s < $ns; ++$s) {
    my $si = $bases->[$s];
    my($all_blank, $pos, $base) = map {$si->{$_}} ('all_blank', 'pos', 'base');
    # don't count indels
    if (!$all_blank && ($pos >= 0) && ($base ne '-')) {
      $snp_counts[$s]++;
    }
  }
}

my $n_snps = scalar(@$snps);
print STDERR "INFO - $n_indels/$n_snps SNP(s) contain at least one indel\n";

for (my $s = 0;$s < $ns; ++$s) {
  my $strain = $strains->[$s];
  my $count = $snp_counts[$s];
}

# order strains by snp count and then alphabetically
my $strain_indices = [0 .. $ns];
my @sorted_strain_indices = sort { ($snp_counts[$b] <=> $snp_counts[$a]) || ($strains->[$a] cmp $strains->[$b]) } @$strain_indices;

for (my $s = 0;$s < $ns; ++$s) {
  my $si = $sorted_strain_indices[$s];
  my $strain = $strains->[$si];
  my $count = $snp_counts[$si];
  print STDERR "INFO - $strain - $count SNP(s)\n";
}

# loop over sorted strains
for (my $s = 0;$s < $ns; ++$s) {
  my $si = $sorted_strain_indices[$s];
  my $strain = $strains->[$si];

  # write VCF file with $strain as reference
  $strain =~ s/\s/_/;
  my $vcf_file = $strain. ".vcf";
  my $vcf_path = File::Spec->catfile($odir, $vcf_file);
  &write_vcf_file($vcf_path, $snp_data, $si);
}

exit(0);

## subroutines
sub check_parameters {
  my $options = shift;
    
  ## make sure required parameters were passed
  my @required = qw(snp_file);
  for my $option ( @required ) {
    unless ( defined $options->{$option} ) {
      die("--$option is a required option");
    }
  }

  ## defaults
  $options->{'output_dir'} = '.' if (!defined($options->{'output_dir'}));

  ## additional parameter checking
  die "SNP file $options->{'snp_file'} not found" if (!-e $options->{'snp_file'});
  my $od = $options->{'output_dir'};
  die "Output directory $od does not exist." if (! -e $od);
  die "Output directory $od is not writable." if (! -w $od);
}

# read SNPs from CloVR-Comparative SNP file
sub read_clovr_snps {
  my($snp_file) = @_;
  my $snp_data = { 'file' => $snp_file };

  # NOTE - the following code is based on Circleator::Parser::SNP_Table
  my $fh = FileHandle->new();
  $fh->open($snp_file)|| die "unable to read from SNP table file $snp_file";
  my $lnum = 0;

  # parse header line
  # e.g., strain1:mol strain1:pos strain1:base strain2:mol strain2:pos strain2:base etc.
  my $header_line = <$fh>;
  ++$lnum;

  # header in example file has an extra column, but this eliminates it
  chomp($header_line);
  my @hf = split(/\t/, $header_line);
  my $n_hcols = scalar(@hf);
  print STDERR "INFO - header line of $snp_file has $n_hcols columns\n";

  # number of columns should always be evenly divisible by 3
  my $num_strains = int($n_hcols / 3);
  if (($num_strains * 3) != $n_hcols) {
    die "header line of $snp_file has $n_hcols columns, which is not evenly divisible by 3";
  }

  # extract strain names from header line
  my $strain_names = [];
  my $col = 0;
  for (my $s = 0;$s < $num_strains;++$s) {
    my $mol_h = $hf[$col++];
    my($strain) = ($mol_h =~ /^(.*):mol$/);
    die "unexpected value at column $col of header line: $mol_h" if (!defined($strain));
    push(@$strain_names, $strain);
    my $pos_h = $hf[$col++];
    die "unexpected value at column $col of header line: $mol_h, not ${strain}:pos" if ($pos_h ne ($strain . ":pos"));
    my $base_h = $hf[$col++];
    die "unexpected value at column $col of header line: $mol_h, not ${strain}:base" if ($base_h ne ($strain . ":base"));
  }

  $snp_data->{'strains'} = $strain_names;
  $snp_data->{'n_strains'} = $num_strains;

  # parse the actual SNP data
  my $snps = [];
  my $num_snps = 0;
  my $num_skipped = 0;
  my $num_snps_with_negative_coords = 0;

 LINE:
  while (my $line = <$fh>) {
    ++$lnum;
    # skip blank lines
    next if ($line =~ /^\s*$/);
    my @f = split(/\t/, $line);
    chomp($f[-1]);
    my $nf = scalar(@f);
    die "wrong number of columns at line $lnum of $snp_file: expected $n_hcols but found $nf" if ($nf != $n_hcols);
    my $col = 0;
    my $snp = [];
    ++$num_snps;
    my $n_negative_coords = 0;
    my $is_indel = 0;

    for (my $s = 0;$s < $num_strains;++$s) {
      my $mol = $f[$col++];
      die "unexpected molecule id '$mol' at column $col of line $lnum" unless ($mol =~ /^(\S.*|)$/);
      my $pos = $f[$col++];
      # negative coordinates showed up in some of the earlier SNP output
      ++$n_negative_coords if ($pos < 0);

      # allow (questionable) negative coordinates
      die "unexpected position '$pos' at column $col of line $lnum" unless ($pos =~ /^(-?\d+|)$/);
      my $base = $f[$col++];
      die "unexpected base '$base' at column $col of line $lnum" unless ($base =~ /^([ACGTUMRWSYKVHDBN\.\-]|)$/i);
      $is_indel = 1 if ($base =~ /\-/);

      # if any 1 of the 3 fields for a strain is blank then they must all be blank
      my $all_blank = 0;
      if (($mol eq '') && ($pos eq '') && ($base eq '')) {
        $all_blank = 1;
      } else {
        die "missing molecule id for $strain_names->[$s] at line $lnum" if ($mol eq '');
        die "missing position for $strain_names->[$s] at line $lnum" if ($pos eq '');
        die "missing base for $strain_names->[$s] at line $lnum" if ($base eq '');
      }
      push(@$snp, { 'all_blank' => $all_blank, 'mol' => $mol, 'pos' => $pos, 'base' => uc($base) });
    }
    push(@$snps, {'is_indel' => $is_indel, 'bases' => $snp });
    ++$num_snps_with_negative_coords if ($n_negative_coords > 0);
  }

  $snp_data->{'snps'} = $snps;
  my $num_parsed = $num_snps - $num_skipped;
  print STDERR "INFO - parsed $num_parsed/$num_snps SNPs from $snp_file\n";
  print STDERR "INFO - $num_snps_with_negative_coords SNPs had negative coordinates\n";

  return $snp_data;
}

sub write_vcf_file {
  my($file, $snp_data, $ref_strain_index) = @_;
  my $fh = FileHandle->new();
  my($ns, $strains, $snps) = map { $snp_data->{$_} } ('n_strains', 'strains', 'snps');
  $fh->open(">$file") || die "unable to write to $file";

  # output VCF meta-information lines
  $fh->print("##fileformat=VCFv4.2\n");
  $fh->print("##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of Samples With Data\">\n");
  $fh->print("##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n");

  # VCF header line
  my @header_cols = @$VCF_COLS;
  # add all strains except reference to the header
  for (my $rsi = 0;$rsi < $ns;++$rsi) {
    next if ($rsi == $ref_strain_index);
    push(@header_cols, $strains->[$rsi]);
  }
  $fh->print(join("\t", @header_cols) . "\n");

  # sort by reference sequence and position - required by VCF spec
  my $rb = sub {
    my $snp = shift;
    my $ref_base = $snp->{'bases'}->[$ref_strain_index];
    return $ref_base;
  };
  my @sorted_snps = sort { my($ba,$bb) = map {&$rb($_)} ($a,$b); ($ba->{'mol'} cmp $bb->{'mol'}) || ($ba->{'pos'} <=> $bb->{'pos'})} @$snps;

  # VCF variants - SNPs
  foreach my $snp (@sorted_snps) {
    my($is_indel, $bases) = map {$snp->{$_}} ('is_indel', 'bases');
    next if ($is_indel);

    my $snp_cols = [];

    # pull out reference base and check whether $all_blank
    my $ref_base = $bases->[$ref_strain_index];
    next if ($ref_base->{'all_blank'});
    my $nb = scalar(@$bases);

    # reference sequence info
    push(@$snp_cols, $ref_base->{'mol'}); # #CHROM
    push(@$snp_cols, $ref_base->{'pos'}); # POS
    push(@$snp_cols, '.'); # ID
    push(@$snp_cols, $ref_base->{'base'}); # REF

    # make allele list
    my $allele_h = {};
    my $anum = 0;
    my $alleles = [];
    # NS - number of samples with data
    my $num_samples = 0;

    for (my $rsi = 0;$rsi < $nb;++$rsi) {
      next if ($rsi == $ref_strain_index);
      my $b = $bases->[$rsi];
      my($all_blank, $mol, $pos, $base) = map {$b->{$_}} ('all_blank', 'mol', 'pos', 'base');
      next if ($all_blank);
      ++$num_samples;
      next if ($base eq $ref_base->{'base'});
      my $index = $allele_h->{$base};
      # add new allele if not present
      if (!defined($index)) {
        $allele_h->{$base} = $index = $anum++;
        push(@$alleles, $base);
      }
      $b->{'allele_index'} = $index;
    }    

    push(@$snp_cols, join(',', @$alleles)); # ALT - comma-delimited list of non-ref alleles using only ACGTN* or <ID>
    push(@$snp_cols, '.'); # QUAL
    push(@$snp_cols, '.'); # FILTER
    push(@$snp_cols, "NS=${num_samples}"); # INFO - semicolon-delimited key=value pair list
    push(@$snp_cols, 'GT'); # FORMAT - colon-delimited list of data types

    for (my $rsi = 0;$rsi < $nb;++$rsi) {
      next if ($rsi == $ref_strain_index);
      my $b = $bases->[$rsi];
      my($all_blank, $mol, $pos, $base, $allele_index) = map {$b->{$_}} ('all_blank', 'mol', 'pos', 'base', 'allele_index');

      # genotype, 0 = reference
      my $gtype = 0;
      if ($all_blank) {
        $gtype = '.';
      }
      elsif ($base ne $ref_base->{'base'}) {
        $gtype = $allele_h->{$base};
        die "FATAL - couldn't find non-reference allele '$base'" if (!defined($gtype));
        # add one because reference is 0 
        $gtype++;
      }
      push(@$snp_cols, $gtype); # genotype - must come first
    }

    # print SNP line
    $fh->print(join("\t", @$snp_cols) . "\n");
  }

  $fh->close();
}
