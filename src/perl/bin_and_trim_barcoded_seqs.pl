#!/usr/bin/perl

=head1 NAME

bin_and_trim_barcoded_seqs.pl - Take a set of barcoded sequences, group them by barcode and trim primers.

=head1 SYNOPSIS

bin_and_trim_barcoded_seqs.pl 
         --input_file=/path/to/input_seqs.fasta
         --reverse_primer='CATGCTGCCTCCCGTAGGAGT'
         --forward_primer='CTGAGCCAGGATCAAACTCT'
         --barcode_file=/path/to/barcodes.txt
        [--output_dir=/path/to/output_dir
         --min_length=200
         --max_edit_dist=2
         --max_barcode_offset=4
         --trim=all
         --log=/path/to/some.log
         --debug=4
         --help ]

=head1 OPTIONS

B<--input_file, -i>
    a multi-FASTA file containing the barcoded sequences to bin and trim.  currently the 
    input sequences MUST be in all-uppercase.

B<--reverse_primer, -r>
    the reverse primer sequence to search for at the beginning of each sequence (following
    the barcode).  if a linker sequence appears between the sequence barcode and the reverse
    primer it should be included in the --reverse_primer sequence.  the script assumes that
    the first two bases of --reverse_primer are a dinucleotide linker sequences, but the 
    only place where this information is used is in a special-case check to see how many 
    sequences contain the barcode and the reverse primer, but NOT the first two characters
    of the reverse primer (the presumed linker sequence)

B<--forward_primer, -f>
    the (reverse complement of) the forward primer sequence to search for at the end of each 
    sequence.

B<--barcode_file,-b>
    a whitespace-delimited text file that specifies the barcodes used to identify sequences in
    --input_file.  each line of the file has two columns: the first specifies the (unique)
    DNA sequence barcode and the second specifies the corresponding well number or sample ID, 
    as needed.  for example:

    ACACACTG    A01
    ACACGTCA    A02
    ACAGACAG    A03
    ACAGCTCA    A04
    
B<--output_dir, -o>
	optional.  directory where the output files will be stored.  this directory will be created
    if it does not exist.  if this option is not specified then the directory that contains 
    --input_file will be used for the output. 
    
B<--min_length, -m>
    optional.  minimum sequence length required after trimming.  (default is 200 bp)

B<--max_edit_dist, -D>
    maximum edit distance between reverse primer and sequence. (default is 2)

B<--max_barcode_offset, -B>
    maximum sequence position at which to look for the (start of) the barcode sequence. (default is 4)
    for example, to ensure that the barcode always starts at the very beginning of each sequence, set
    --max_barcode_offset=0

B<--trim, -t>
    optional. specifies what to trim off the sequences:
      'all' (the default) - trim both the barcode and the primer(s)
      'barcodes' - trim only the barcode but leave the primer(s) intact
      'none' - bin the sequences based on barcode but do not do any barcode or primer trimming
    Note that the forward primer may not always be present (e.g., when the read length of the 
    underlying sequencing technology is significantly shorter than the expected insert size.)

B<--log,-l>
    optional.  path to a log file the script should create.  will be overwritten if
    it already exists.

B<--debug,-d>
    optional.  the debug level for the logger (an integer)

B<--help,-h>
    This help message/documentation.

=head1  DESCRIPTION

Take a set of barcoded sequences, group them by barcode and trim primers.  The script
takes as input a single multi-FASTA file of sequences and divides them amongst a set of
output files based on the persence of an embedded sequence barcode (e.g., "ACACACTG").
The script will trim the barcode sequence off each sequence, in addition to trimming 
the reverse and forward primers (--reverse_primer and --forward_primer), if they can
be found.  The --trim option can be used to specify the trimming behavior: trim neither
the barcode nor the primers, trim only the barcode, or trim both the primers and the
barcode.

Each sequence with a recognizable barcode will be placed in an output file whose name
contains both the barcode and also the well/sample ID specified for that barcode in 
the --barcode_file.  For example, if the input sequence file is "input.fsa", the barcode
detected in the sequence is "TGCATCGA", and the sample id that corresponds to that 
barcode is "H12" then the (trimmed) sequence will be placed in one of the following 
two files:

input.fsa.H12.TGCATCGA.filtered
input.fsa.H12.TGCATCGA.discarded

A sequence will be placed in the ".filtered" file if both the barcode and the reverse
primer can be found (within --max_edit_dist), and the length of the sequence after 
trimming is greater than or equal to the minimum value specified by --min_length.  
If the reverse primer cannot be found or the length of the sequence after trimming is
below the threshold then the (perhaps partially trimmed) sequence will be placed in 
the ".discarded" file.  (If reverse primer trimming is disabled then the presence of
the reverse primer isn't required for the sequence to be placed in the ".filtered"
file.)  Any sequence without a recognizable barcode will be placed in a discard
file named based solely on the input file, e.g.:

test.fsa.discarded

The FASTA header line of the sequences (in both the .filtered and .discarded files) 
will be updated to reflect the trimming and barcode/primer matching.  For example, 
in the following defline:

>000123_456_789 length=261 uaccno=XXFABC40343|E01|CTGAGTGT|trimmed_length:232|rev_primer_mismatches:0|fwd_primer_mismatches:NA

The new fields are:

E01: 
 well/sample id that corresponds to the barcode 'CTGAGTGT'
CTGAGTGT: 
 literal barcode that was identified in the sequence (and trimmed, unless --trim=none)
trimmed_length
 the length of the sequence after trimming the barcode and/or primers
rev_primer_mismatches: 
 edit distance between the canonical reverse primer given by --reverse_primer and the 
 primer found in the sequence.  should always be <= --max_edit_dist.  "NA" means that
 the primer was not found.
fwd_primer_mismatches:
 edit distance between the canonical forward primer given by --forward_primer and the 
 primer found in the sequence.  this value will always be either 0 (exact match) or 
 "NA" (no match) since approximate matching of forward primers is not yet supported.

=head1  INPUT

A multi-FASTA formatted file of sequences to bin and trim (specified by --input_file) 
and a set of unique barcodes (specified by --barcode_file).  Currently the input sequences,
barcodes, and primers must be all-uppercase.

=head1  OUTPUT

A set of multi-FASTA formatted files as described above.  The script also prints the 
following information to STDOUT:

Levenshtein edit distance histogram for reverse RNA primer:

    distance      num_seqs
           0        239058
           1          5014
           2*          289

       total        244361
* - setting for --max_edit_dist

Indel count histogram for reverse RNA primer:

      indels      num_seqs
          -2           155
          -1          1925
           0        241723
           1           503
           2            55

       total        244361

Number of reverse primers that match exactly except for missing CA linker: 45

=head1  CONTACT

    William Hsiao
    william.hsiao@gmail.com

=cut

use strict;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Pod::Usage;
use Text::LevenshteinXS qw(distance);
use File::Basename;
use File::Spec;
use Ergatis::Logger;

## global/default values
my $DEFAULT_MIN_LENGTH = 200;
my $DEFAULT_MAX_EDIT_DIST = 2;
my $DEFAULT_MAX_BARCODE_OFFSET = 4;
my $DEFAULT_TRIM = 'all';

my $MAX_BARCODE_LENGTH = undef;

# edit distance histogram for reverse primer
my $REV_EDIT_DIST_HIST = {};

# indel histogram for reverse primer e.g.,
# -2: 2 deletions in input seq relative to $reverse_primer
# +2: 2 insertions in input seq relative to $reverse_primer
my $REV_INDEL_HIST = {};

# count of sequences missing the dinucleotide linker
my $NUM_MISSING_DINUC       = 0;

# overall sequence counts
my $NUM_SEQS = 0;
my $NUM_SEQS_WITH_BARCODE = 0;
my $NUM_SEQS_WITH_BARCODE_AND_REV_PRIMER = 0;

# cache of open filehandles for writing binned/trimmed sequences
my $WRITE_FH_CACHE = {};

## input/options
my $options = {};
my $results = GetOptions($options, 
                         'input_file|i=s', 
                         'barcode_file|b=s',
                         'reverse_primer|r=s',
                         'forward_primer|f=s',
                         'output_dir|o=s',    
                         'min_length|m=i',
                         'max_edit_dist|D=i', 
                         'max_barcode_offset|B=i', 
                         'trim|t=s',
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

## MAIN SECTION
my $barcodes = get_barcode_plate_mapping($logger, $options->{'barcode_file'});
my $input_file = $options->{'input_file'};
my $input_basefilename = basename($input_file);
my $output_dir = $options->{'output_dir'};
$output_dir =~ s/\/$//;
my $output_file_prefix = File::Spec->catfile($output_dir, $input_basefilename);

if (!-e $output_dir) {
    my $res = mkdir $output_dir;
    die "failed to mkdir $output_dir: $!";
}

#warn user and quit if any output files already exist
check_for_existing_output_files($output_file_prefix, $barcodes);

my $ifh;
if ( $input_file =~ /\.gz$/ ) {
    open( $ifh, "<:gzip", $input_file )
        || die "can't read file $input_file: $!";
}
else {
    open( $ifh, "<$input_file" ) || die "can't read file $input_file: $!";
}

# get the sequences from the file
my $seq = '';
my @headers;
my $is_first = 1;

while (<$ifh>) {
    chomp;
    if (/^\>/) {
        push @headers, $_;
        if ($is_first) {
            $is_first = 0;
        }
        else {
            check_and_reset( $options, $logger, $barcodes, $output_file_prefix, $headers[-2], \$seq );
        }
    }
    else {
        $seq .= $_;
    }
}

# check the last one
check_and_reset( $options, $logger, $barcodes, $output_file_prefix, $headers[-1], \$seq );

# print edit distance histogram
print "\nLevenshtein edit distance histogram for reverse primer:\n\n";
my @dists = sort { $a <=> $b } keys %$REV_EDIT_DIST_HIST;
if (scalar(@dists) == 0) { 
    print " ERROR - reverse primer not found in any sequence\n";
} else {
    &print_histogram( $REV_EDIT_DIST_HIST, 0, $dists[-1], 'distance',
                      'num_seqs', { $options->{max_edit_dist} => '*' }, 
                      $NUM_SEQS_WITH_BARCODE_AND_REV_PRIMER );
    print "* - setting for --max_edit_dist\n";
}

# print indel count histogram
print "\nIndel count histogram for reverse primer:\n\n";
my @nindels = sort { $a <=> $b } keys %$REV_INDEL_HIST;
if (scalar(@nindels) == 0) { 
    print " ERROR - reverse primer not found in any sequence\n";
} else {
    &print_histogram( $REV_INDEL_HIST, $nindels[0], $nindels[-1], 'indels',
                      'num_seqs', {}, $NUM_SEQS_WITH_BARCODE_AND_REV_PRIMER );
    print "\nNumber of reverse primers that match exactly except for missing dinucleotide linker: $NUM_MISSING_DINUC\n";
}
print "\n";

&close_write_fhs();
exit(0);

## subroutines

# this is used to avoid doing open/close for every sequence
#
# $barcodes - barcode hashref
# $prefix - file path prefix used to name output files
# $barcode - a valid key in $barcodes
# $type - either 'filtered' or 'discarded'
#
sub get_write_fh {
    my($barcodes, $prefix, $barcode, $type) = @_;
    my $filename = defined($barcode) ? 
        "$prefix.$barcodes->{$barcode}.$barcode.$type" : "$prefix.$type";
    my $fh = $WRITE_FH_CACHE->{$filename};
    if (!defined($fh)) {
        die "$filename already exists" if (-e $filename);
        open($fh, ">$filename") || die "unable to write to $filename: $!";
        $WRITE_FH_CACHE->{$filename} = $fh;
    }
    return $fh;
}

# close all filehandles previously opened by get_write_fh
#
sub close_write_fhs {
    foreach my $file (keys %$WRITE_FH_CACHE) {
        my $fh = $WRITE_FH_CACHE->{$file};
        if (close $fh) {
            delete $WRITE_FH_CACHE->{$file};
        } else {
            die "error closing filehandle for $file";
        }
    }
}

# process a single FASTA sequence in the input file
#
# $options - command-line options
# $logger - Ergatis::Logger
# $barcodes - barcode hashref
# $prefix - file path prefix used to name output files
# $header - FASTA header line, including '>'
# $seq_ref - reference to the nucleotide sequence
#
sub check_and_reset {
	my ( $options, $logger, $barcodes, $prefix, $header, $seq_ref ) = @_;
    # remove whitespace
	$$seq_ref =~ s/\s//g;
    # error if sequence is not uppercase (since matching is done on uppercase primers only)
	die "the following sequence is not all-uppercase: $header"
	  if ( $$seq_ref ne uc($$seq_ref) );
    ++$NUM_SEQS;

    my $orig_leng = length($$seq_ref);
	my $barcode = check_for_barcode($options, $logger, $barcodes, $header, $seq_ref);

	if ( defined($barcode) ) {
		my $trimseq  = undef;
		my $seq_leng = 0;
        my $rev_edit_dist = undef;
        my $fwd_edit_dist = undef;
        ++$NUM_SEQS_WITH_BARCODE;

		# don't trim barcode or primers if --trim=none
		if ( $options->{'trim'} =~ /none/i ) {
			$trimseq  = $$seq_ref;
			$seq_leng = length($trimseq);
		}
        # trim barcode, track # bases trimmed
		else {
			$trimseq = trim_barcode( $logger, $header, $seq_ref, $barcode, $options->{'max_barcode_offset'} );
            # this should never fail since check_for_barcode has already confirmed the barcode's presence
            die "internal error: check_barcode succeeded, but trim_barcode failed" if (!defined($trimseq));
			$seq_leng = length($trimseq);

			# trim reverse and forward primers unless --trim=barcodes|none
			if ( $options->{'trim'} !~ /barcodes/ ) {
				($trimseq, $rev_edit_dist) = trim_reverse_primer( $options, $logger, $header, \$trimseq, uc($options->{'reverse_primer'}) );
                ++$NUM_SEQS_WITH_BARCODE_AND_REV_PRIMER if (defined($trimseq));
				$seq_leng = length($trimseq);

                # trim forward primer, if present
                my ($f_trim, $f_ed) = trim_forward_primer( $header, \$trimseq, uc($options->{'forward_primer'}));
                if (defined($f_trim)) {
                    $trimseq = $f_trim;
                    $fwd_edit_dist = $f_ed;
                    $seq_leng = length($trimseq);
                }
            }
        }

		my $above_lengthcutoff = ( $seq_leng >= $options->{'min_length'} );
        my $rpm = "rev_primer_mismatches:" . (defined($rev_edit_dist) ? $rev_edit_dist : 'NA');
        my $fpm = "fwd_primer_mismatches:" . (defined($fwd_edit_dist) ? $fwd_edit_dist : 'NA');
        my $new_header = join('|', $header, $barcodes->{$barcode}, $barcode, "trimmed_length:${seq_leng}", $rpm, $fpm);

		if ( ( defined($trimseq) ) && $above_lengthcutoff ) {
            my $ofh = &get_write_fh($barcodes, $prefix, $barcode, 'filtered');
            print $ofh "${new_header}\n${trimseq}\n";
		}
		elsif ( !defined($trimseq) ) {
            my $ofh = &get_write_fh($barcodes, $prefix, $barcode, 'discarded');
            print $ofh "${new_header}\n$${seq_ref}\n";
		}
		else {
            my $ofh = &get_write_fh($barcodes, $prefix, $barcode, 'discarded');
            print $ofh "${new_header}\n${trimseq}\n";
		}
	}
	else {
        $logger->info("sequence '$header' does not have a unique barcode: adding it to $prefix.discarded");
        my $nobc_discardfh = &get_write_fh($barcodes, $prefix, undef, 'discarded');
		print $nobc_discardfh "$header\n$$seq_ref\n";
	}
	$$seq_ref = '';
}

# determine which barcode(s) a sequence matches
#
# $options - command-line options
# $logger - Ergatis::Logger
# $barcodes - barcode hashref
# $header - FASTA header line, including '>'
# $seq_ref - reference to the nucleotide sequence
#
sub check_for_barcode {
	my ($options, $logger, $barcodes, $header, $seq_ref) = @_;

    # the barcode must appear within this range for it to count as a match:
    my $max_barcode_end_offset = $options->{'max_barcode_offset'} + $MAX_BARCODE_LENGTH;
	my $primer_seq = substr( $$seq_ref, 0, $max_barcode_end_offset );
    my $matches = {};

    $logger->debug("checking for barcode in seq '$primer_seq'");

	foreach my $bc ( keys %$barcodes ) {
		if ( $primer_seq =~ /$bc/g ) {
            $matches->{$bc} = pos($primer_seq) - length($bc);
		}
	}
    
    my @matching_barcodes = keys %$matches;
    my $nm = scalar(@matching_barcodes);

    if ($nm == 1) {
        return $matching_barcodes[0];
    } elsif ($nm == 0) {
        return undef;
    } else {
        # TODO - pick the first match here
        my @sorted = sort { $matches->{$a} <=> $matches->{$b} } keys %$matches;
        my $ms = join(',', map { $_ . "(position=" . $matches->{$_} . ")" } @sorted);
        $logger->warn("sequence '$header' matched $nm barcodes in first $options->{'max_barcode_offset'} bp, choosing the first: $ms");
        return $sorted[0];
    }
}

# check that none of the output files exists; fail if so
#
# $prefix - file path prefix used to name output files
# $barcodes - barcode hashref
#
sub check_for_existing_output_files {
	my($prefix, $barcodes) = @_;
    
    # global discard file
	if ( -e "$prefix.discarded" ) {
		die "$prefix.discarded already exists\n";
	}

    # barcode-specific files
	foreach my $barcode ( keys %$barcodes ) {
        my @bfiles = map { 
            my $file = join('.', $prefix, $barcodes->{$barcode}, $barcode, $_);
            die "$file already exists" if (-e $file);
        } 
        ('filtered', 'discarded');
	}
}

# Take a sequence and a barcode and trim off the barcode and any sequences before it
# If the barcode is not found at or before $max_offset a warning is printed (and no
# trimming is done)
#
# $logger - Ergatis::Logger
# $header - FASTA header line, including '>'
# $seq_ref - reference to the nucleotide sequence
# $barcode - barcode to trim from $$seq_ref
# $max_offset - maximum (starting) sequence offset at which to match the barcode
#
sub trim_barcode {
	my ( $logger, $header, $seq_ref, $barcode, $max_offset ) = @_;
	my $trimmed_seq = undef;
	$$seq_ref =~ /$barcode/g;
	$$seq_ref =~ /\G(\w+)/;
    my $match_posn = pos($$seq_ref);

    if (defined($match_posn)) {
        $trimmed_seq = $1;
        my $max_end_offset = $max_offset + length($barcode);
    
        if ($match_posn > $max_end_offset) {
            my $bop = $match_posn - length($barcode);
            $logger->warn("sequence '$header' matches the barcode '$barcode' at position $bop, which is greater than the --max_barcode_offset");
            $trimmed_seq = undef;
        }
    }
	return $trimmed_seq;
}

# trim --forward_primer from the sequence.  the subroutine expects to find the forward 
# primer towards the end of the sequence (or not at all) and will trim the primer and 
# everything that follows it
#
# $header - FASTA header line, including '>'
# $seq_ref - reference to the nucleotide sequence 
# $primer - forward primer sequence to trim
#
sub trim_forward_primer {
	my ( $header, $seq_ref, $primer ) = @_;
    my $trimmed_seq = undef;
    my $edit_dist = undef;

    # check for exact match
    if ($$seq_ref =~ /^(.*)$primer/) {
        $trimmed_seq = $1;
        $edit_dist = 0;
    }

    # TODO - check for inexact matches
    return ($trimmed_seq, $edit_dist);
}

# trim --reverse_primer from the sequence.  the subroutine expects that the barcode has
# already been identified and trimmed.
#
# $options - command-line options
# $logger - Ergatis::Logger
# $header - FASTA header line, including '>'
# $seq_ref - reference to the nucleotide sequence
# $primer - reverse primer sequence to trim
#
sub trim_reverse_primer {
	my ( $options, $logger, $header, $seq_ref, $primer ) = @_;
	my $trimmed_seq = undef;
    my $edit_dist = undef;
	my $primer_seq = substr( $$seq_ref, 0, length($primer) );

	# check for exact match
	if ( $primer eq $primer_seq ) {
        $edit_dist = 0;
		$REV_EDIT_DIST_HIST->{0}++;
		$REV_INDEL_HIST->{0}++;
		$trimmed_seq = substr( $$seq_ref, length($primer) );
        return ($trimmed_seq, $edit_dist);
	} 

	# special case - check for missing dinucleotide linker 
    # (the first two bases of the --reverse_primer)
	my $primer_nd = $primer;
	$primer_nd =~ s/^\S\S//i;
	my $seq_primer_nd = substr( $$seq_ref, 0, length($primer) - 2 );
	if (( $options->{'max_edit_dist'} >= 2 ) && ( $primer_nd eq $seq_primer_nd )) {
		$REV_EDIT_DIST_HIST->{2}++;
		$REV_INDEL_HIST->{-2}++;
		$trimmed_seq = substr( $$seq_ref, length($primer_nd) );
		$NUM_MISSING_DINUC++;
        $edit_dist = 2;
	}

    # if not an exact match then try searching for the best nearby match
    # i.e., do a local search around length($primer) and pick the minimum edit distance
	else {
		my $radius = $options->{'max_edit_dist'};
		my $dists  = [];
		my $rpl    = length($primer);

		for ( my $r = -$radius ; $r <= $radius ; ++$r ) {
			my $sslen = $rpl + $r;
			$primer_seq = substr( $$seq_ref, 0, $sslen );

			# compute Levenshtein edit distance
			my $ldist = distance( $primer, $primer_seq );
			push(
				@$dists,
				{
					'dist'   => $ldist,
					'abs_r'  => abs($r),
					'r'      => $r,
					'len'    => $sslen,
					'primer' => $primer_seq
				}
			);
		}

		# err on the side of overtrimming by a base or two
		my @sortedDists = sort {
			     ( $a->{'dist'} <=> $b->{'dist'} )
			  || ( $a->{'abs_r'} <=> $b->{'abs_r'} )
			  || ( $b->{'len'} <=> $a->{'len'} )
		} @$dists;

        $logger->debug("sorted Levenshtein edit distances for sequence '$header':");
        foreach my $sd (@sortedDists) {
            my $dd = Data::Dumper->new( [ $sd ] );
            $dd->Indent(0);
            $logger->debug(" dist=" . $dd->Dump());
        }

		my $best = $sortedDists[0];
		if ( $best->{'dist'} <= $options->{'max_edit_dist'} ) {
			$trimmed_seq = substr( $$seq_ref, $best->{'len'} );
            $edit_dist = $best->{'dist'};
			$REV_EDIT_DIST_HIST->{ $best->{'dist'} }++;
			$REV_INDEL_HIST->{ $best->{'r'} }++;
			my $primer = $best->{'primer'};

            $logger->debug( "reverse primer match summary:");
            $logger->debug( "   header: $header" );
            $logger->debug( "     dist: " . $best->{'dist'});
            $logger->debug( "        r: " . $best->{'r'} );
            $logger->debug( "      seq: $$seq_ref" );
            $logger->debug( "   primer: $primer" );
            $logger->debug( "   primer: $primer" );

			$primer =~ s/\S/ /g;
            $logger->debug( "  trimmed: ${primer}${trimmed_seq}");
		}
		else {
            $logger->info( "reverse primer sequence not found in sequence '$header', within max_edit_dist=" . $options->{'max_edit_dist'});
		}
	}
	return ($trimmed_seq, $edit_dist);
}

# get mapping from barcode to well/sample id
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
        $logger->logdie("barcode file contains multiple rows for '$id'") if (defined($barcodes{uc($id)}));
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
    if ($nb == 0) {
        $logger->logdie("no barcodes read from $barcode_file");
    } else {
        $logger->info("$nb barcode(s) read from $barcode_file");
    }
    $logger->error("max_barcode_length = 0") if ($MAX_BARCODE_LENGTH == 0);
    $logger->warn("max_barcode_length < 6") if ($MAX_BARCODE_LENGTH < 6);

	return \%barcodes;
}

sub check_parameters {
	my $options = shift;

    ## make sure required parameters were passed
    my @required = qw(input_file reverse_primer forward_primer barcode_file);
    for my $option ( @required ) {
        unless  ( defined $options->{$option} ) {
            die "--$option is a required option";
        }
    }

	## defaults
	$options->{'min_length'} = $DEFAULT_MIN_LENGTH if (!defined($options->{'min_length'}));
	$options->{'max_edit_dist'} = $DEFAULT_MAX_EDIT_DIST if (!defined($options->{'max_edit_dist'}));
	$options->{'max_barcode_offset'} = $DEFAULT_MAX_BARCODE_OFFSET if (!defined($options->{'max_barcode_offset'}));
	$options->{'output_dir'} = dirname($options->{'input_file'}) if (!defined($options->{'output_dir'}));
    $options->{'trim'} = $DEFAULT_TRIM if (!defined($options->{'trim'}));

    my $to = $options->{'trim'};
    die "illegal --trim option (${to})" unless ($to =~ /^all|barcodes|none$/i);

    my $if =  $options->{'input_file'};
    if ((!-e $if) || (!-r $if)) {
        die "$if does not exist or is not readable";
    }

    my $bf =  $options->{'barcode_file'};
    if ((!-e $bf) || (!-r $bf)) {
        die "$bf does not exist or is not readable";
    }

    if ($options->{'max_barcode_offset'} < 0) {
        die "illegal value for --max_barcode_offset";
    }
}

# print a simple histogram of values
#
# $hist - the histogram, encoded as a hashref
# $min - minimum integer key value to print, inclusive
# $max - maximum integer key value to print, inclusive
# $keyName - name for the hashref keys
# $valName - name for the hashref values
# $highlightValues - hashref mapping key to highlight to a single character (e.g., '*', '!')
# $expectedSum - check sum against this value (if defined) and complain if it doesn't match
#
sub print_histogram {
	my ( $hist, $min, $max, $keyName, $valName, $highlightValues, $expectedSum ) = @_;
    my $sum = 0;
	printf( "%12s%1s %12s\n", $keyName, '', $valName );
	for ( my $i = $min ; $i <= $max ; ++$i ) {
		my $val = $hist->{$i};
		$val = 0 if ( !defined($val) );
		my $hc = $highlightValues->{$i} || '';
		printf( "%12s%1s %12s\n", $i, $hc, $val );
        $sum += $val;
	}
    printf( "\n%12s%1s %12s\n", 'total', '', $sum );
    if (defined($expectedSum) && ($sum != $expectedSum)) {
        $logger->error("histogram sum ($sum) does not match expected value ($expectedSum)\n");
    }
}

=head2 TODO:
-Check additional barcodes beyond the first (in the case of multiple matches) if choosing
 the first one fails to match the primer sequence(s)
-Relax requirement that all the input sequences and primers have to be all-uppercase.
-Allow approximate matches to the forward primer.
-Check the product size in cases where both primers are found?
