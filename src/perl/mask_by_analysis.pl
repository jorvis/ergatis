#!/usr/local/bin/perl

use lib (@INC,$ENV{"PERL_MOD_DIR"});
no lib "$ENV{PERL_MOD_DIR}/i686-linux";
no lib ".";

=head1 NAME
mask_by_analysis.pl  - mask FASTA sequences using results from an analysis

=head1 SYNOPSIS

USAGE:  mask_by_analysis.pl -i analysis_bsml.bsml -o masked_sequence_bsml.bsml -a analysisname,... 
                         [-f featureclass,...] [-x N] [-l /path/to/logfile.log]

=head1 OPTIONS

B<--input,-i>
    The input BSML file with sequences to be masked and analysis to mask with.

B<--output,-o>
    The output path.

B<--output_bsml,-b>
	Output BSML file name.

B<--mask_char,-x>
    The character to mask with (default = X).

B<--analysis_types,-a>
    Comma or whitespace delimitted list of Analysis types for masking (eg: wu-blastp).
	
B<--feature_types,-f>
    Comma or whitespace delimitted list of Feature classes for masking (eg: exon).
	[optional]

B<--random>
    Mask NA sequences with random bases (0 = no, 1 = yes: default 0).
	
B<--log,-l>
    Log file

B<--help,-h>
    This help message

=head1   DESCRIPTION

This script can be used to mask a fasta sequence file based on the results of
some analysis. All Seq-pair-alignment regions, and/or Feature interval-loc
regions produced by a given analysis method can be masked on the input sequence.

NOTE:

Calling the script name with no options or --help will display the syntax requirements.

=head1  OUTPUT

The sequence produces fasta format sequence files for each sequence on which
masking is performed. An accompanying BSML document is produced which describes
these sequence files.

=head1  CONTACT
        Brett Whitty
        bwhitty@tigr.org

=cut

use strict;
use warnings;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use Pod::Usage;
use XML::Twig;
use Ergatis::Logger;
use BSML::BsmlBuilder;

$| = 1;
my %options = ();
GetOptions (\%options,
            'input|i=s',
			'analysis_types|a=s',
			'feature_types|f:s',
			'mask_char|x=s',
			'output|o=s',
			'output_bsml|b=s',
			'random:i',
			'log|l=s',
			'debug|d=i',
            'help|h') || pod2usage();


my $infile;
my $outfile;
my $mask_char;
my $mask_regions;
my $bsml_sequences;
my $analysis_types; 
my $feature_types;
my $global_counter = 0;

if ($options{'help'}){
	pod2usage(verbose => 2);
}

my $logfile = $options{'log'} || Ergatis::Logger::get_default_logfilename();
my $logger = new Ergatis::Logger('LOG_FILE'=>$logfile,
                                  'LOG_LEVEL'=>$options{'debug'});
$logger = $logger->get_logger();

if (!$options{'input'}) {
	pod2usage("bsml input file name must be provided with --input");
}
if (!-e $options{'input'}) {
	if (-e $options{'input'}.".gz") {
		$infile = $options{'input'}.".gz";
	} else {
		$logger->logdie("input file '$options{input}' doesn't exist");
	}
} else {
	$infile = $options{'input'};
}

if (!$options{'output'}) {
	pod2usage("must specify output path with --output");
}
if (!$options{'output_bsml'}) {
	pod2usage("must specify output bsml file name with --output_bsml");
}

$options{'output'} =~ s/\/$//;
unless (-d $options{'output'}) {
	pod2usage("argument to --output '$options{output}' is not a directory");
}

if (!$options{'mask_char'}) {
	$mask_char = 'X';
} else {
	if (length($options{'mask_char'}) > 1) {
		pod2usage("mask character '$options{mask_char}' is more than one character");
	}
	$mask_char = $options{'mask_char'};
}

if (!$options{'analysis_types'}){
	pod2usage("at least one analysis type to mask on must be provided with --analysis_types");
}
if ($options{'analysis_types'}) {
	$analysis_types = do_flag_hash($options{'analysis_types'});
}
if ($options{'feature_types'}) {
	$feature_types = do_flag_hash($options{'feature_types'});
}

$outfile = $options{'output'}."/".$options{'output_bsml'};

my $infh;
my $outfh;

srand();

## find all interval-locs and seq-pair-runs
my $loc_twig = XML::Twig->new(
							twig_handlers => {
								'Seq-pair-alignment' => \&seqpair_handler,
								'Feature'			 => \&feature_handler,
											 }
						);

if ($infile =~ /\.gz$/) {
	open($infh, "<:gzip", $infile) 
	 || logger->logdie("couldn't open gzipped input file '$infile'");
} else {
   	open($infh, $infile)
	 || logger->logdie("couldn't open input file '$infile'");
}	
$loc_twig->parse($infh);
close $infh;

## mask the sequence an build a new BSML document
my $seq_twig = new XML::Twig(	TwigRoots => {
		'Sequence' => 1,
											 },
								TwigHandlers => {'Sequence' => \&sequence_handler}
							);
					 
$seq_twig->set_pretty_print('indented');
 
my $doc = new BSML::BsmlBuilder();

if ($infile =~ /\.gz$/) {
	open($infh, "<:gzip", $infile) 
	 || logger->logdie("couldn't open gzipped input file '$infile'");
} else {
   	open($infh, $infile)
	 || logger->logdie("couldn't open input file '$infile'");
}	
$seq_twig->parse($infh);
close $infh;

$doc->write($outfile);

exit();

## deals with seq-pair-alignments
## pulls out coordinates for seq-pair-runs
## and stores them in mask_regions
sub seqpair_handler {
	my ($twig, $seq_pair_alignment) = @_;

	my $seq_id = $seq_pair_alignment->{'att'}->{'refseq'};
	
	my $flag = 0;
	
	## choose whether we want to deal with this sequence element
	## based on the contents of its analysis link
	my @links = $seq_pair_alignment->children('Link');
	foreach my $link(@links) {
		if ($link->{'att'}->{'rel'} eq 'analysis') {
			## may want to make role a flag option
			if ($link->{'att'}->{'role'} eq 'computed_by') {
				$link->{'att'}->{'href'} =~ /^#(.*)_analysis$/;
				my $analysis_id = $1;
				if ($analysis_types->{$analysis_id}) {
					$flag = 1;
				}
			}
		}
	}

	unless ($flag) { 
		return;
	}
	
	foreach my $seq_pair_run($seq_pair_alignment->children('Seq-pair-run')) {
		push(@{$mask_regions->{$seq_id}}, 
			 [$seq_pair_run->{'att'}->{'refpos'}, 
			  $seq_pair_run->{'att'}->{'runlength'}]
		    ); 
	}
}

## deals with feature elements
## pulls out the interval-locs
## and stores them in mask_regions
sub feature_handler {
	my ($twig, $feat) = @_;
    my $id = $feat->att('id');
	
	my $parent_seq = $feat->parent('Sequence');
	my $seq_id = $parent_seq->{'att'}->{'id'};
	
	my $flag = 0;
	
	## choose whether we want to deal with this sequence element
	## based on the contents of its analysis link
	my @links = $feat->children('Link');
	foreach my $link(@links) {
		if ($link->{'att'}->{'rel'} eq 'analysis') {
			## may want to make role a flag option
			if ($link->{'att'}->{'role'} eq 'computed_by') {
				$link->{'att'}->{'href'} =~ /^#(.*)_analysis$/;
				my $analysis_id = $1;
				if ($analysis_types->{$analysis_id}) {
					$flag = 1;
				}
			}
		}
	}

	unless ($flag) { 
		return;
	}

	## only process specified feature classes if any were specified
	if ($options{'feature_types'} && !$feature_types->{$feat->att('class')}) {
		return;
	}
	
    my $seq_int = $feat->first_child('Interval-loc');
		
    my ($start_pos, $end_pos) = ($seq_int->att('startpos'), $seq_int->att('endpos'));
		
	if ($start_pos > $end_pos) {
		($end_pos, $start_pos) = ($start_pos, $end_pos);
	}
		
	my $len = $end_pos - $start_pos;
		
	push(@{$mask_regions->{$seq_id}}, 
		 [$start_pos, 
		  $len]
	    );
}
		 

## deals with sequence elements and makes sure that
## we have sequence for the sequences we want to mask
## then it does the masking and writes modified BSML
## to the output file
sub sequence_handler {
	my ($twig, $sequence) = @_;
	
	my $seq_id;
	my $class;
	my $molecule;
	my $seq_data_import;
	my $seq_data;
	my $source;
	my $identifier;
	my $format;
	
	
	my $flag = 0;
	
	## choose whether we want to deal with this sequence element
	## based on the contents of its analysis link
	my @links = $sequence->children('Link');
	foreach my $link(@links) {
		if ($link->{'att'}->{'rel'} eq 'analysis') {
			## may want to make role a flag option
			if ($link->{'att'}->{'role'} eq 'input_of') {
				if ($link->{'att'}->{'href'} =~ /^#(.*)_analysis$/) {
					my $analysis_id = $1;
					if ($analysis_types->{$analysis_id}) {
						$flag = 1;
					}
				} else {
					$logger->logdie("couldn't parse analysis type");
				}
			}
		}
	}

	unless ($flag) { 
		return;
	}

	$seq_id = $sequence->{'att'}->{'id'};
	$class = $sequence->{'att'}->{'class'};
	$molecule = $sequence->{'att'}->{'molecule'};
	
	my $sequence_file = $options{'output'}."/$seq_id.mask_by_analysis.fsa";
	
	if (!defined($class)) {
		if ($logger->is_debug()) {$logger->debug("WARNING: sequence class of '$seq_id' was not defined\n");}
	} else {
		if ($options{'random'} && $class eq 'polypeptide') {
			$logger->logdie("Can't use option --random with polypeptide sequences");
		}
	}
	if (!defined($molecule) || $molecule eq 'mol-not-set') {
		if ($logger->is_debug()) {$logger->debug("WARNING: molecule type of '$seq_id' was not defined\n");}
	} else {
		if ($options{'random'} && $molecule eq 'aa') {
			$logger->logdie("Can't use option --random with polypeptide sequences");
		}
	}
	
	if ($seq_data_import = $sequence->first_child('Seq-data-import')) {
		## sequence is referenced via seq-data-import

		$source = $seq_data_import->{'att'}->{'source'};
		$identifier = $seq_data_import->{'att'}->{'identifier'};
		$format = $seq_data_import->{'att'}->{'format'};
		
		## we only support pulling sequences from fasta files right now
		## but get_sequence_by_id could be modified to support other formats 
		if (defined($format) && $format ne 'fasta') {
			$logger->logdie("unsupported seq-data-import format '$format' found");
		}
		
		unless (-e $source) {
			$logger->logdie("fasta file referenced in BSML Seq-data-import '$source' doesn't exist");
		}
		unless (defined($identifier)) {
			$logger->logdie("Seq-data-import for '$seq_id' does not have a value for identifier");
		}
		
		my $seq = get_sequence_by_id($source, $identifier);
			
		if (length($seq) > 0) {

			$seq = mask_sequence($seq, $mask_regions->{$seq_id}, $mask_char);
			
			write_seq_to_fasta($sequence_file, $seq_id, $seq);

			$seq_data_import->{'att'}->{'source'} = $sequence_file;
			
		} else {
			$logger->logdie("couldn't fetch sequence for '$seq_id' from Seq-data-import source '$source'");
		}
		
	} elsif ($seq_data = $sequence->first_child('Seq-data')) {
		## sequence is in the BSML
		my $seq = mask_sequence($seq_data->text(), $mask_regions->{$seq_id}, $mask_char);
		
		write_seq_to_fasta($sequence_file, $seq_id, $seq);
		
		$seq_data->delete();
		
		my $seq_data_import = new XML::Twig::Elt('Seq-data-import', '');
		$seq_data_import->set_atts({
										'format' 		=>	'fasta',
										'id'			=> 	'_Bsml_'.$global_counter++,
										'identifier'	=> 	$seq_id,
										'source'		=> 	$sequence_file,
							 	   });
    	$seq_data_import->paste( 'last_child', $sequence);         
		
	} else {
		## there is no Seq-data or Seq-data-import for the sequence
		$logger->logdie("No sequence present in BSML sequence element");
	}

	my $seq_imp = $sequence->first_child('Seq-data-import');
	
	my $bsml_seq = $doc->createAndAddExtendedSequenceN(
						%{$sequence->{'att'}},
												      );
	$doc->createAndAddSeqDataImportN(
					'seq' => $bsml_seq,
					%{$seq_imp->{'att'}},
									);
	$doc->createAndAddLink(
							$bsml_seq, 
							'analysis',
							'#mask_by_analysis_analysis',
							'computed_by',
						  );
}

## pull a sequence from a fasta file by sequence id
## where the sequence id is the header string up to
## the first whitespace char
sub get_sequence_by_id {
	my ($fname, $id) = @_;
	my $seq_id = '';
	my $sequence = '';
	open (IN, $fname) || $logger->logdie("couldn't open fasta file for reading");
	TOP: while (<IN>) {
		chomp;
		if (/^>([^\s]+)/) {
			$seq_id = $1;
			if ($seq_id eq $id) {
				while (<IN>) {
					chomp;
					if (/^>/) {
						last TOP;
					} else {
						$sequence .= $_;
					}
				}
			}	
		}
	}
	close IN;

	return $sequence;
}

## writes a nicely formatted fasta file
sub write_seq_to_fasta {
	my ($file, $header, $sequence) = @_;
	
	open (OUT, ">$file") || $logger->logdie("couldn't write fasta file '$file'");
		
	$sequence =~ s/\W+//g;
	$sequence =~ s/(.{1,60})/$1\n/g;
		
	print OUT ">$header\n$sequence";
	close OUT;
}

## does the dirty business of masking the sequence
sub mask_sequence {
	my ($seq, $segments, $mask_char) = @_;

	## not gonna do any interval merging here because I'm lazy
	## and hopefully these string operations should be fast
	foreach my $segment(@{$segments}) {
		my ($start, $len) = @{$segment};
		if ($start < 0) {
			$logger->logdie("mask_sequence was passed start < 0: start=$start len=$len\n");
		}
		if ($len < 0) {
			$logger->logdie("mask_sequence was passed len < 0: start=$start len=$len\n");
		}
		if (($start + $len) > length($seq)) {
			$logger->logdie("mask_sequence was passed start + len > length(seq): start=$start len=$len seqlen=".length($seq)."\n");
		}
		
		## if the random flag, mask with random NTs
		if ($options{'random'}) {
			$seq = substr($seq, 0, $start).random_nt($len).substr($seq, $start + $len);
		## else we'll mask with mask_char
		} else {
			$seq = substr($seq, 0, $start).($mask_char x $len).substr($seq, $start + $len);
		}
	}
	
	return $seq;
}

## populates a hash with flags for a set of keys
sub do_flag_hash{
	my ($i) = @_;
	my $hash;
	
	my @vals = split(/[,;:\s]+/, $i);
	foreach my $val(@vals) {
		$hash->{$val} = 1;
	}
	
	return $hash;
}

## returns a string of random nucleotides of determinate length
sub random_nt {
	my ($len) = @_;
	my $seq = '';
	my $nt = 'ATGC';
	
	for (my $i=1; $i<= $len; $i++) {
		my $val = int(rand(3)) + 1;
		$nt = substr($nt, $val) . substr($nt, 0, $val);
		$seq .= substr($nt, int(rand(4)), 1);
	}
	return $seq;
}
