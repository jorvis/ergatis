#!/usr/local/bin/perl
BEGIN{foreach (@INC) {s/\/usr\/local\/packages/\/local\/platform/}};
use lib (@INC,$ENV{"PERL_MOD_DIR"});
no lib "$ENV{PERL_MOD_DIR}/i686-linux";
no lib ".";

=head1  NAME 

iprscan2bsml.pl - convert iprscan btab output to BSML

=head1 SYNOPSIS

USAGE: iprscan2bsml.pl 
    --input=/path/to/somefile.iprscan.raw 
    --output=/path/to/somefile.iprscan.bsml
  [ --fasta_input=/path/to/iprscan/input.fsa
    --gzip_output=1 ]

=head1 OPTIONS

B<--input,-i> 
    Input tab-delimited file from an iprscan search.

B<--query_id,-r>
	ID of the query sequence.  This allows the creation
    of a Sequence element in BSML even if there are no matches.

B<--query_file_path,-q>
	Full path to FASTA file containing query sequence.  This allows the creation
    of a Sequence element in BSML even if there are no matches.

B<--debug,-d> 
    Debug level.  Use a large number to turn on verbose debugging. 

B<--gzip_output,-g>
    Optional.  Non-zero value will give you a compressed bsml file. If output
    does not contain a .gz extension, it will be added automatically.

B<--log,-l> 
    Log file

B<--output,-o> 
    Output BSML file

B<--help,-h> 
    This help message

=head1   DESCRIPTION

This script is used to convert the tab-delimited output from a nap search into BSML.

=head1 INPUT

The input to this component is the btab file created by the iprscan executable.  This file
is 14-column, tab-delimited text document.  The input must be a nap btab output from
a single search.  Fromthe documentation, the raw output looks like this:

    ------
    NF00181542      0A5FDCE74AB7C3AD        272     HMMPIR  PIRSF001424     Prephenate dehydratase  1       270     6.5e-141        T       06-Oct-2004         IPR008237       Prephenate dehydratase with ACT region  Molecular Function:prephenate dehydratase activity (GO:0004664), Biological Process:L-phenylalanine biosynthesis (GO:0009094)
    ------

	Where: NF00181542:             is the id of the input sequence.
	       27A9BBAC0587AB84:       is the crc64 (checksum) of the proteic sequence (supposed to be unique).
	       272:                    is the length of the sequence (in AA).
	       HMMPIR:                 is the anaysis method launched.
	       PIRSF001424:            is the database members entry for this match.
	       Prephenate dehydratase: is the database member description for the entry.
	       1:                      is the start of the domain match.
	       270:                    is the end of the domain match.
	       6.5e-141:               is the evalue of the match (reported by member database anayling method).
	       T:                      is the status of the match (T: true, ?: unknown).
	       06-Oct-2004:            is the date of the run.
	       IPR008237:              is the corresponding InterPro entry (if iprlookup requested by the user).
	       Prephenate dehydratase with ACT region:                           is the description of the InterPro entry.
	       Molecular Function:prephenate dehydratase activity (GO:0004664):  is the GO (gene ontology) description for the InterPro entry.


Illegal characters will be removed from the 
IDs for the query sequence (column 0) and subject hit (column 5) if necessary to create 
legal XML id names.  For each element, the original, unmodified name will be stored in the "title"
attribute of the Sequence element.  You should make sure that your ids don't begin with
a number.  This script will successfully create a BSML document regardless of your ID names,
but the resulting document may not pass DTD validation.

=head1 OUTPUT

The BSML file to be created is defined using the --output option.  If the file already exists
it will be overwritten.

=head1 CONTACT

    Joshua Orvis
    jorvis@tigr.org

=cut

use strict;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use Pod::Usage;
use Ergatis::Logger;
use BSML::BsmlRepository;
use BSML::BsmlBuilder;
use BSML::BsmlParserTwig;
use File::Basename;
use XML::Twig;

my $gzip;
my $fasta_input;
my $defline;
my $identifier;

my %options = ();
my $results = GetOptions (\%options, 
			  'input|i=s',
              'output|o=s',
              'fasta_input|f=s',
              'gzip_output|g=s',
			  'query_id|r=s',
			  'query_file_path|q=s',
              'log|l=s',
              'debug=s',
			  'help|h') || pod2usage();

my $logfile = $options{'log'} || Ergatis::Logger::get_default_logfilename();
my $logger = new Ergatis::Logger('LOG_FILE'=>$logfile,
				  'LOG_LEVEL'=>$options{'debug'});
$logger = $logger->get_logger();

# display documentation
if( $options{'help'} ){
    pod2usage( {-exitval=>0, -verbose => 2, -output => \*STDOUT} );
}

## make sure all passed options are peachy
&check_parameters(\%options);

## we want a new doc
my $doc = new BSML::BsmlBuilder();

## open the input file for parsing
open (my $ifh, $options{'input'}) || $logger->logdie("can't open input file for reading");

my %seqs_found;

## add the query sequence to the BSML doc
my $seq = $doc->createAndAddSequence($options{query_id}, $options{query_id}, undef, 'na', 'assembly');
$doc->createAndAddSeqDataImport($seq, 'fasta', $options{query_file_path}, '', $options{query_id});
$seq->addBsmlLink('analysis', '#iprscan_analysis', 'input_of');
$doc->createAndAddBsmlAttribute( $seq, 'defline', $defline) if($defline);
$seqs_found{$options{query_id}} = 1;

while (<$ifh>) {
    ## ignore whitespace lines
    next if ( /^\s*$/ );
    chomp;
    
    ## there should be 14 elements in cols, unless we have an unrecognized format.
    my @cols = split("\t");
    unless (scalar @cols >= 13) {
        $logger->error("the following iprscan line was not recognized and could not be parsed (should have 14 columns, actually has " . scalar(@cols) . "):\n$_\n") if ($logger->is_error);
        next;
    }
    
    ## This isn't the right way to do the subject ID, but given a name like PIRSF001424 there
    ##  isn't a reliable way to parse the database name from the entry name, unless we have a 
    ##  fixed list of database prefixes we can explicitly check.  Just do this for now so we
    ##  can have a working solution.
    my ($qry_id, $sbj_db, $sbj_id, $sbj_name) = ($cols[0], $cols[3], $cols[4], $cols[5]);
    
    ## the qry ID only counts up to the first whitespace
    if ($qry_id =~ /(.+?)\s+/) {
        $qry_id = $1;
    }
    
    ## make sure both of these are valid IDs
    $qry_id =~ s/[^a-zA-Z0-9\.\-\_]/_/g;
    $sbj_id =~ s/[^a-zA-Z0-9\.\-\_]/_/g;
    
    ## has this subject sequence been added to the doc yet?
    if (! exists $seqs_found{$sbj_id}) {
        my $seq = $doc->createAndAddSequence($sbj_id, $cols[5], undef, 'aa', 'polypeptide');
        $doc->createAndAddCrossReference( parent => $seq, 
                                          database => $sbj_db, 
                                         'identifier-type' => 'accession',
                                          identifier => $sbj_id );
        $doc->createAndAddCrossReference( parent => $seq, 
                                          database => $sbj_db, 
                                         'identifier-type' => 'entry-name',
                                          identifier => $sbj_name );
        if ($cols[11] =~ /IPR/ ) {
            $doc->createAndAddCrossReference( parent => $seq, 
                                              database => 'interpro', 
                                             'identifier-type' => 'accession',
                                              identifier => $cols[11] );
        }
        $seqs_found{$sbj_id} = 1;
    }

    ## each line represents a matched domain, so create a Seq-pair-alignment
    ## skipped attributes here are complength, compstart, compend and method
    my $seq_pair_aln = $doc->createAndAddSequencePairAlignment( refseq => $qry_id,
                                                                refstart => $cols[6] - 1,
                                                                refend => $cols[7],
                                                                reflength => $cols[2],
                                                                compseq => $sbj_id,
                                                                class => 'match',
                                                              );
    $seq_pair_aln->addBsmlLink('analysis', '#iprscan_analysis');

    if ( $cols[8] ne 'NA' ) {
        ## add the e_value (will be the same for each matching segment)
        $doc->createAndAddBsmlAttribute($seq_pair_aln, 'e_value', $cols[8]);
    }
    
    ## now add the Seq-pair-run
    ## skipped attributes are runprob
    my $run = $doc->createAndAddSequencePairRun(   alignment_pair => $seq_pair_aln,
                                                   runscore => 0,  ## raw score is not given
                                                   runlength => $cols[7] - $cols[6],
                                                   comprunlength => $cols[7] - $cols[6],
                                                   refpos => $cols[6] - 1,
                                                   refcomplement => 0,
                                                   comppos => 0,
                                                   compcomplement => 0,
                                               );

    $doc->createAndAddBsmlAttribute($run, 'class', 'match_part');
    
    if ( $cols[8] ne 'NA' ) {
        $doc->createAndAddBsmlAttribute($run, 'e_value', $cols[8]);
    }

}

## add the analysis element
my $analysis = $doc->createAndAddAnalysis(
                                            id => 'iprscan_analysis',
                                            sourcename => $options{'output'},
                                         );
$doc->createAndAddBsmlAttribute( $analysis, 'version', 'current' );
$doc->createAndAddBsmlAttribute( $analysis, 'algorithm', 'iprscan' );

## now write the doc
$doc->write($options{'output'}, '', $gzip);

exit;

sub check_parameters {
    
    ## make sure input file exists
    if (! -e $options{'input'}) { $logger->logdie("input file $options{'input'} does not exist") }

    ## make user an output file was passed
    if (! $options{'output'}) { $logger->logdie("output option required!") }
    
    unless ($options{query_file_path} && $options{query_id}) { 
        $logger->logdie("query_file_path and query_id options required!") 
    }
    
    ## do we want compressed output?
    if($options{'gzip_output'}) {
        $gzip = 1;
    } else {
        $gzip = 0;
    }    

    ## store some options if the fasta_input file was included
    if($options{'query_file_path'}) {
        $fasta_input = $options{'query_file_path'};
        open(IN, "< $fasta_input") or
            $logger->logdie("Could not open $fasta_input ($!)");
        while(<IN>) {
            chomp;
            if(/^>(.*)/) {
                $defline = $1;
                $identifier = $1 if($defline =~ /^([^\s]+)/);
            }
        }
        close(IN);
    } 

    return 1;
}

