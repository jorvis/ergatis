#!/usr/bin/perl

=head1 NAME

count_barcoded_seqs.pl - Summarize the output of bin_and_trim_barcoded_seqs.pl in tabular form.

=head1 SYNOPSIS

count_barcoded_seqs.pl
         --input_dir=/path/to/bin_and_trim_barcoded_seqs/output_dir
         --input_file=/path/to/input_seqs.fasta
         --barcode_file=/path/to/barcodes.txt
         --output_file=/path/to/summary-report.txt
        [--log=/path/to/some.log
         --debug=4
         --help ]

=head1 OPTIONS

B<--input_dir, -D>
    directory into which bin_and_trim_barcoded_seqs.pl wrote the output sequence files

B<--input_file,-i>
    path to the same --input_file used by bin_and_trim_barcoded_seqs.pl

B<--barcode_file,-b>
    path to the same --barcode_file used by bin_and_trim_barcoded_seqs.pl

B<--output_file,-o>
    path to the file in which to write the tabular summary

B<--log,-l>
    optional.  path to a log file the script should create.  will be overwritten if
    it already exists.

B<--debug,-d>
    optional.  the debug level for the logger (an integer)

B<--help,-h>
    This help message/documentation.

=head1  DESCRIPTION

Summarize the output of bin_and_trim_barcoded_seqs.pl in tabular form.

=head1  INPUT

The output of bin_and_trim_barcoded_seqs.pl.

=head1  OUTPUT

A table listing, among other things, the number of trimmed and discarded sequences for each barcode.

=head1  CONTACT

    Jonathan Crabtree
    jonathancrabtree@gmail.com

=cut

use strict;
use FileHandle;
use File::Spec;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Pod::Usage;
use Ergatis::Logger;

## globals
my $MAX_BARCODE_LENGTH = undef;

## input/options
my $options = {};
my $results = GetOptions($options, 
                         'input_dir|D=s', 
                         'input_file|i=s', 
                         'barcode_file|b=s', 
                         'output_file|o=s', 
                         'debug|d=s',
                         'log|l=s',
                         'help|h'
                         ) || pod2usage();

## display documentation
if ( $options->{'help'} ) {
	pod2usage( { -exitval => 0, -verbose => 2, -output => \*STDERR } );
}

&check_parameters($options);

## initialize logging
my $logfile = $options->{'log'} || Ergatis::Logger::get_default_logfilename();
my $logger = new Ergatis::Logger('LOG_FILE'=>$logfile, 'LOG_LEVEL'=>$options->{'debug'});
$logger = Ergatis::Logger::get_logger();

## main program
my $barcodes = &get_barcode_plate_mapping($logger, $options->{'barcode_file'});
my $nb = scalar(keys %$barcodes);
$logger->info("read $nb barcode(s) from $options->{'barcode_file'}");

# the alternate output file suffixes (".out" and ".discard") are used by older 
# versions of the binning/trimming script:
opendir(DDIR, $options->{'input_dir'});
my @files = grep(/\.(out|discard|filtered|discarded)$/, readdir(DDIR));
closedir(DDIR);

# determine shared file prefix
my $prefix = undef;
foreach my $file (@files) {
    if ((!defined($prefix) && ($file =~ /^(.*)\.[^\.]+\.[^\.]+\.(out|discard|filtered|discarded)$/))) {
        $prefix = $1;
        last;
    }
}

die "unable to determine file prefix" if (!defined($prefix));
$logger->debug("file prefix=$prefix");

# sequence counts
my $total_num_seqs = 0;
my $total_num_out = 0;
my $total_num_discard = 0;
# counts of seqs with no reverse primer detected
my $total_num_out_nrp = 0;
my $total_num_discard_nrp = 0;
# counts of seqs with no forward primer detected
my $total_num_out_nfp = 0;
my $total_num_discard_nfp = 0;

# count number of sequences in input
my ($num_in, $num_in_nrp, $num_in_nfp) = &count_multi_FASTA_seqs($options->{'input_file'});

# output file
my $ofh = FileHandle->new();
my $of = $options->{'output_file'};
$ofh->open(">$of") || die "unable to write to $of";

# print table header
printf $ofh "%12s %12s %12s %12s %16s %16s %16s %16s\n", 
    'BARCODE', 'NUM_OUT', 'NUM_DISCARD', 'DISCARD_%', 'NO_REV_PRIMER', 
    'NO_REV_PRIMER_%', 'NO_FWD_PRIMER', 'NO_FWD_PRIMER_%';

# barcode-specific counts
foreach my $barcode (sort keys %$barcodes) {
    my $sample = $barcodes->{$barcode};
    my $out = File::Spec->catfile($options->{'input_dir'}, join('.', $prefix, $sample, $barcode, 'filtered'));
    my $discard = File::Spec->catfile($options->{'input_dir'}, join('.', $prefix, $sample, $barcode, 'discarded'));

    # check for older file suffixes:
    if ((!-e $out) && (!-e $discard)) {
        $out = File::Spec->catfile($options->{'input_dir'}, join('.', $prefix, $sample, $barcode, 'out'));
        $discard = File::Spec->catfile($options->{'input_dir'}, join('.', $prefix, $sample, $barcode, 'discard'));
    }

    $logger->debug("reading from output file $out");
    my ($num_out, $num_out_nrp, $num_out_nfp) = -e $out ? &count_multi_FASTA_seqs($out) : (0,0,0);
    $logger->debug("reading from discard file $discard");
    my ($num_discard, $num_discard_nrp, $num_discard_nfp) = -e $discard ? &count_multi_FASTA_seqs($discard) : (0,0,0);

    # any sequence in the regular output file should match the reverse primer
    $logger->warn("sequences with no reverse primer found in $out") if ($num_out_nrp > 0);

    $total_num_out += $num_out;
    $total_num_discard += $num_discard;
    $total_num_out_nrp += $num_out_nrp;
    $total_num_discard_nrp += $num_discard_nrp;
    $total_num_out_nfp += $num_out_nfp;
    $total_num_discard_nfp += $num_discard_nfp;
    $total_num_seqs += ($num_out + $num_discard);

    # any output file should have at least one sequence
    $logger->warn("file $out exists but has no sequences") if (-e $out && ($num_out == 0));
    $logger->warn("file $discard exists but has no sequences") if (-e $discard && ($num_discard == 0));
    
    my $sum = $num_out + $num_discard;
    my $discard_pct = ($sum == 0) ? '0' : sprintf("%0.1f", ($num_discard/$sum) * 100.0);
    my $num_nrp = $num_discard_nrp + $num_out_nrp;
    my $num_nfp = $num_discard_nfp + $num_out_nfp;
    my $nrp_pct = ($sum == 0) ? '0' : sprintf("%0.1f", ($num_nrp/$sum) * 100.0);
    my $nfp_pct = ($sum == 0) ? '0' : sprintf("%0.1f", ($num_nfp/$sum) * 100.0);

    printf $ofh  "%12s %12s %12s %12s %16s %16s %16s %16s\n", 
           $barcode, $num_out, $num_discard, $discard_pct . '%', $num_nrp, $nrp_pct . '%', $num_nfp, $nfp_pct . '%';
}

# count sequences with no recognizable barcode (those in the general discard file)
my $no_barcode = File::Spec->catfile($options->{'input_dir'}, join('.', $prefix, 'discarded'));

# also try ".discard" for reverse compatibility
if (!-e $no_barcode) {
    $no_barcode = File::Spec->catfile($options->{'input_dir'}, join('.', $prefix, 'discard'));
}

my ($num_no_barcode, $num_no_rp, $num_no_fp) = &count_multi_FASTA_seqs($no_barcode);
$total_num_seqs += $num_no_barcode;
print $ofh "\n";

my $total_num_nrp = $total_num_discard_nrp + $total_num_out_nrp;
my $total_num_nfp = $total_num_discard_nfp + $total_num_out_nfp;
printf $ofh  "%12s %12s %12s %12s %16s %16s %16s\n", 'SUBTOTALS:', 
    $total_num_out, $total_num_discard, '', $total_num_nrp, '', $total_num_nfp;
printf $ofh "%12s %12s %12s %12s %16s\n", 'NO BARCODE:', '0', $num_no_barcode, '', '0';

print $ofh "\n total input seqs: $num_in\n";
print $ofh "total output seqs: $total_num_seqs\n";

my $exitval = 0;
if ($num_in != $total_num_seqs) {
    print $ofh "ERROR: sequence count mismatch!";
    $logger->error("total number of sequences in $options->{'input_file'} ($num_in) does not match total number of output sequences ($total_num_seqs)");
    $exitval = 1;
}

$ofh->close();
exit($exitval);

## subroutines

sub check_parameters {
	my $options = shift;

    ## make sure required parameters were passed
    my @required = qw(input_dir input_file barcode_file output_file);
    for my $option ( @required ) {
        unless  ( defined $options->{$option} ) {
            die "--$option is a required option";
        }
    }

    my $id =  $options->{'input_dir'};
    if ((!-r $id) || (!-d $id)) {
        die "$id is not readable or is not a directory";
    }

    my $if =  $options->{'input_file'};
    if ((!-e $if) || (!-r $if)) {
        die "$if does not exist or is not readable";
    }

    my $bf =  $options->{'barcode_file'};
    if ((!-e $bf) || (!-r $bf)) {
        die "$bf does not exist or is not readable";
    }
}

# get mapping from barcode to well/sample id
# (copied from bin_and_trim_barcoded_seqs.pl)
#
# $logger - Ergatis::Logger
# $barcode_file - full path to file that contains barcode - well/sample ID mapping
#
sub get_barcode_plate_mapping {
    my($logger, $barcode_file) = @_;
	my %barcodes;
    my $different_lengths = 0;

    # parse barcode file
    open( my $if, $barcode_file )
        or die "Unable to open the barcode file - $barcode_file";
    while (<$if>) {
        chomp;
        my ( $id, $well ) = split( /\s+/, $_ );
        # keys _must_ be in all-caps as all matching is done in uppercase
        $barcodes{uc($id)} = $well;
        my $bl = length($id);
        $MAX_BARCODE_LENGTH = $bl if (!defined($MAX_BARCODE_LENGTH) || ($bl > $MAX_BARCODE_LENGTH));
        if (defined($MAX_BARCODE_LENGTH) && ($MAX_BARCODE_LENGTH != $bl)) {
            $different_lengths = 1;
        }
    }
    close $if;
    $logger->warn("barcode lengths are not uniform") if ($different_lengths);
    my $nb = scalar(keys %barcodes);
    $logger->warn("no barcodes read from $barcode_file") if ($nb == 0);
    $logger->error("max_barcode_length = 0") if ($MAX_BARCODE_LENGTH == 0);
    $logger->warn("max_barcode_length < 6") if ($MAX_BARCODE_LENGTH < 6);

	return \%barcodes;
}

# count the number of sequences in a multi-FASTA file
# 
# $file - path to a multi-FASTA file
#
sub count_multi_FASTA_seqs {
    my($file) = @_;
    my $nseqs = 0;
    my $nnrp = 0; # num with no reverse primer
    my $nnfp = 0; # num with no forward primer
    my $fh = FileHandle->new();
    $fh->open($file, 'r') || die "unable to read from $file";
    while (my $line = <$fh>) {
        if ($line =~ /^>/) {
            ++$nseqs;
            # old defline format
            ++$nnrp if ($line =~ /no primer/i);
            # new defline format
            ++$nnrp if ($line =~ /rev_primer_mismatches:\s*NA/i);
            ++$nnfp if ($line =~ /fwd_primer_mismatches:\s*NA/i);
        }
    }
    $fh->close();

    return ($nseqs, $nnrp, $nnfp);
}
