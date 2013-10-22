#!/usr/bin/perl

=head1  NAME 

remove-nonunique-clovr-comparative-snps.pl - Remove all SNPs from a CloVR-Comparative SNP file _except_ those unique to a specified reference.

=head1 SYNOPSIS

remove-nonunique-clovr-comparative-snps.pl
         --snp_file=somefile_mini.snps
         --reference=Bordetella_holmesii_41130_BH003
        [--help
         --man]

=head1 OPTIONS

B<--snp_file,-s>
    Path to a SNP output file prodcued by an instance of the CloVR-Comparative ergatis pipeline.

B<--reference,-r>
    The reference genome/taxon: only SNPs unique to this genome will be retained in the output.

B<--help,-h> 
    Display the script documentation.

B<--man,-m>
    Display the script documentation.

=head1 DESCRIPTION

Remove all SNPs from a CloVR-Comparative SNP file _except_ those
unique to a specified reference.  Used by the CloVR-Comparative
pipeline to reduce the size of the input SNP file before parsing it
and converting it to the space-inefficient BioPerl representation.

=head1 INPUT

The SNP output file prodcued by an instance of the CloVR-Comparative ergatis pipeline.

=head1 OUTPUT

The same SNP file but with all SNPs not unique to the given reference taxon removed.

=head1 CONTACT

    Jonathan Crabtree
    jonathancrabtree@gmail.com

=cut

use strict;
use FileHandle;
use Getopt::Long;

## input
my $options = {};

&GetOptions($options, 
            "snp_file|s=s",
            "reference|r=s",
            "help|h",
            "man|m"
           );

## display documentation
if ( $options->{'help'} || $options->{'man'} ) {
  pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDERR} );
}
&check_parameters($options);

## input
my($snp_file, $ref_strain) = map { $options->{$_} } ('snp_file', 'reference');

## main program
my $fh = FileHandle->new();
$fh->open($snp_file)|| die "unable to read from SNP table file $snp_file";
my $lnum = 0;

# parse and echo the header line
# e.g., strain1:mol strain1:pos strain1:base strain2:mol strain2:pos strain2:base etc.
my $header_line = <$fh>;
++$lnum;
print $header_line;

# the header in our example file had an extra column, but this eliminates it
chomp($header_line);
my @hf = split(/\t/, $header_line);
my $n_hcols = scalar(@hf);
print STDERR "header line of $snp_file has $n_hcols columns\n";

# number of columns should be evenly divisible by 3
my $num_strains = int($n_hcols / 3);
if (($num_strains * 3) != $n_hcols) {
  print STDERR "header line of $snp_file has $n_hcols columns, which is not evenly divisible by 3\n";
}

my $strain_names = [];
# an index into $strain_names for the reference strain
my $ref_strain_index = undef;
my $col = 0;
for (my $s = 0;$s < $num_strains;++$s) {
  my $mol_h = $hf[$col++];
  my($strain) = ($mol_h =~ /^(.*):mol$/);
  die("unexpected value at column $col of header line: $mol_h") if (!defined($strain));
  push(@$strain_names, $strain);
  my $pos_h = $hf[$col++];
  die("unexpected value at column $col of header line: $mol_h, not ${strain}:pos") if ($pos_h ne ($strain . ":pos"));
  my $base_h = $hf[$col++];
  die("unexpected value at column $col of header line: $mol_h, not ${strain}:base") if ($base_h ne ($strain . ":base"));
  $ref_strain_index = $s if ($strain eq $ref_strain);
}
die("snp-ref '$ref_strain' not found in $snp_file") if (!defined($ref_strain_index));

# iterate over the rest of the file, filtering out non-unique SNPs
LINE:
while (my $line = <$fh>) {
  ++$lnum;
  # skip blank lines
  next if ($line =~ /^\s*$/);
  my @f = split(/\t/, $line);
  chomp($f[-1]);
  my $nf = scalar(@f);
  die("wrong number of columns at line $lnum of $snp_file: expected $n_hcols but found $nf") if ($nf != $n_hcols);
  my $col = 0;
  my $snp_info = [];
  
  for (my $s = 0;$s < $num_strains;++$s) {
    my $mol = $f[$col++];
    die("unexpected molecule id '$mol' at column $col of line $lnum") unless ($mol =~ /^(\S.*|)$/);
    my $pos = $f[$col++];
    # TODO - find out whether these represent a bug in the underlying SNP pipeline
    if ($pos < 0) {
      print STDERR "negative position '$pos' at column $col of line $lnum, skipping SNP\n";
      next LINE;
    }
    die("unexpected position '$pos' at column $col of line $lnum") unless ($pos =~ /^(\d+|)$/);
    my $base = $f[$col++];
    die("unexpected base '$base' at column $col of line $lnum") unless ($base =~ /^([ACGTUMRWSYKVHDBN\.\-]|)$/i);
    
    # if any 1  of the 3 fields is blank then they must all be blank
    my $all_blank = 0;
    if (($mol eq '') && ($pos eq '') && ($base eq '')) {
      $all_blank = 1;
    } else {
      die("missing molecule id for $strain_names->[$s] at line $lnum") if ($mol eq '');
      die("missing position for $strain_names->[$s] at line $lnum") if ($pos eq '');
      die("missing base for $strain_names->[$s] at line $lnum") if ($base eq '');
    }
    push(@$snp_info, { 'all_blank' => $all_blank, 'mol' => $mol, 'pos' => $pos, 'base' => $base });
  }
  
  # skip if SNP doesn't appear in reference
  my $ref_snp_info = $snp_info->[$ref_strain_index];
  next if ($ref_snp_info->{'all_blank'});
  my($ref_mol, $ref_pos, $ref_base) = map {$ref_snp_info->{$_}} ('mol', 'pos', 'base');

  # count number of non-reference positions that are the same as the reference
  my $num_same_as_ref = 0;
  for (my $s = 0;$s < $num_strains;++$s) {
    next if ($s == $ref_strain_index);
    my $snp_info = $snp_info->[$s];
    my($mol, $pos, $base, $all_blank) = map {$snp_info->{$_}} ('mol', 'pos', 'base', 'all_blank');
    next if ($all_blank);
    ++$num_same_as_ref if (uc($base) eq uc($ref_base));
  }

  # print only unique SNPs
  #
  # Note that this counts a SNP as unique (to the reference) so long
  # as it differs from every other strain in which there is a hit
  # (versus requiring there to be a hit in each non-reference strain
  # AND requiring that those hits are different.) Deciding between
  # these two options is left to the downstream code (i.e., Circleator
  # config file.)
  #
  print $line if ($num_same_as_ref == 0); 
}

exit(0);

## subroutines

sub check_parameters {
  my $options = shift;
    
  ## make sure required parameters were passed
  my @required = qw(snp_file reference);
  for my $option ( @required ) {
    unless ( defined $options->{$option} ) {
      die("--$option is a required option");
    }
  }
  $options->{'reference'} =~ s/(\.gene_summary)$//;
}
