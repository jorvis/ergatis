#!/usr/bin/perl

BEGIN{foreach (@INC) {s/\/usr\/local\/packages/\/local\/platform/}};
use lib (@INC,$ENV{"PERL_MOD_DIR"});
no lib "$ENV{PERL_MOD_DIR}/i686-linux";
no lib ".";

=head1  NAME 

nap2bsml.pl - convert nap btab output to BSML

=head1 SYNOPSIS

USAGE: nap2bsml.pl 
    --input=/path/to/somefile.nap.btab 
    --output=/path/to/somefile.nap.bsml

=head1 OPTIONS

B<--input,-i> 
    Input btab file from a nap search.

B<--query_file_path,-q>
    Full path to FASTA file containing query sequence.

B<--query_id>
    ID of query sequence

B<--cutoff_identity> 
    Filter results on % identity (exclude < cutoff). 

B<--cutoff_similarity> 
    Filter results on % similarity (exclude < cutoff). 
     
B<--debug,-d> 
    Debug level.  Use a large number to turn on verbose debugging. 

B<--log,-l> 
    Log file

B<--output,-o> 
    Output BSML file

B<--qzip_output,-g>
    Optional.  A non-zero value will make compressed output.

B<--help,-h> 
    This help message

=head1   DESCRIPTION

This script is used to convert the btab output from a nap search into BSML.

=head1 INPUT

The input to this component is the btab file created by the nap executable.  This file
is 19-column, tab-delimited text document.  The input must be a nap btab output from
a single search.  Seq-pair-alignment elements are based on the chain ID from nap, so if you concatenate
btabs together each chainID 1, for example, would get merged into the same alignment, which 
is probably not what you want.  Define the input file using the --input option.

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

my $qfDefline;

my %options = ();
my $results = GetOptions (\%options, 
              'input|i=s',
              'query_file_path|q=s',
              'query_id=s',
              'output|o=s',
              'gzip_output|g=s',
              'cutoff_identity:s',
              'cutoff_similarity:s',
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
my $ifh;
if($options{'input'} =~ /\.gz$/) {
    open ($ifh, "<:gzip", $options{'input'}) || $logger->logdie("can't open input file for reading");
} else {
    open ($ifh, $options{'input'}) || $logger->logdie("can't open input file for reading");
}

## each chainID = one Seq-pair-alignment
## each chain segment = one Seq-pair-run

my %seqs_found;

## this hash will hold the chainIDs as their key and a reference to its Seq-pair-alignment
my %chains;
while (<$ifh>) {
    ## ignore whitespace lines
    next if ( /^\s*$/ );
    chomp;
    
    ## there should be 19 elements in cols, unless we have an unrecognized format.
    my @cols = split("\t");
    unless (scalar @cols == 19) {
        $logger->error("the following nap btab line was not recognized and could not be parsed (should have 19 columns, actually has " . scalar(@cols) . "):\n$_\n") if ($logger->is_error);
        next;
    }
       
    if ($options{'cutoff_identity'} && $cols[10] < $options{'cutoff_identity'}) {
        next;
    }
    if ($options{'cutoff_similarity'} && $cols[11] < $options{'cutoff_similarity'}) {
        next;
    }
 
    my ($qry_id, $sbj_id) = ($cols[0], $cols[5]);
    
    ## the qry ID only counts up to the first whitespace
    if ($qry_id =~ /(.+?)\s+/) {
        $qry_id = $1;
    }
    
    ## has this query sequence been added to the doc yet?
    if (! exists $seqs_found{$qry_id}) {
        my $seq = $doc->createAndAddSequence($qry_id, $cols[0], undef, 'na', 'assembly');
        $doc->createAndAddSeqDataImport($seq, 'fasta', $options{'query_file_path'}, '', $cols[0]);
        $seq->addBsmlLink('analysis', '#aat_aa_analysis', 'input_of');
        $seq->addBsmlAttr('defline', $qfDefline) if($qfDefline);
        $seqs_found{$qry_id} = 1;
    }
    
    ## has this subject sequence been added to the doc yet?
    if (! exists $seqs_found{$sbj_id}) {
        my $seq = $doc->createAndAddSequence($sbj_id, $cols[5], undef, 'aa', 'polypeptide');
        $seq->addBsmlAttr( 'defline', "$cols[5] $cols[15]" );
        $doc->createAndAddSeqDataImport($seq, 'fasta', $cols[4], '', $cols[5]);
        $doc->createAndAddCrossReferencesByParse( sequence => $seq, string => $cols[5]);
        $seqs_found{$sbj_id} = 1;
    }
    
    my ($chainID, $segmentID) = ($cols[13], $cols[14]);
    
    ## if this combination doesn't exist yet, create its Seq-pair-alignment
    if (! exists $chains{$chainID}) {
        ## skipped attributes here are complength, compstart, compend, compxref, refxref and method
        $chains{$chainID} = $doc->createAndAddSequencePairAlignment( refseq => $qry_id,
                                                                     refxref => ":$qry_id",
                                                                     refstart => 0,
                                                                     refend => $cols[2] - 1,
                                                                     reflength => $cols[2],
                                                                     compseq => $sbj_id,
                                                                     compxref => "$cols[4]:$sbj_id",
                                                                     class => 'match',
                                                                   );
        $chains{$chainID}->addBsmlLink('analysis', '#aat_aa_analysis', 'computed_by');

        ## add the total_score (will be the same for each matching segment)
        $doc->createAndAddBsmlAttribute($chains{$chainID}, 'total_score', $cols[18]);
    }
    
    ## now add the Seq-pair-run
    ## skipped attributes are runprob
    my $run = $doc->createAndAddSequencePairRun(   alignment_pair => $chains{$chainID},
                                                   runscore => $cols[12],
                                                   runlength => abs($cols[7] - $cols[6]) + 1,
                                                   comprunlength => abs($cols[9] - $cols[8]) + 1,
                                                   refpos => min($cols[6], $cols[7]) - 1,
                                                   refcomplement => $cols[17] eq 'Minus' ? 1 : 0,
                                                   comppos => min($cols[8], $cols[9]) - 1,
                                                   compcomplement => 0,
                                               );
    $doc->createAndAddBsmlAttribute($run, 'class', 'match_part');
    $doc->createAndAddBsmlAttribute($run, 'percent_identity', $cols[10]);
    $doc->createAndAddBsmlAttribute($run, 'percent_similarity', $cols[11]);
    $doc->createAndAddBsmlAttribute($run, 'chain_number', $chainID);
    $doc->createAndAddBsmlAttribute($run, 'segment_number', $segmentID);
}

## if there were no results this will create a sequence stub
my $align = &createAndAddNullResult(
                                        doc             => $doc,
                                        query_name      => $options{'query_id'},
                                        query_length    => '',
                                        class           => 'assembly',
                                   );
     
## add the analysis element
my $analysis = $doc->createAndAddAnalysis(
                            id => 'aat_aa_analysis',
                            sourcename => $options{'output'},
                          );

## now write the doc
$doc->write($options{'output'}, '', $options{'gzip_output'});

exit;


sub check_parameters {
    
    ## make sure input file exists
    if (! -e $options{'input'}) { 
        $logger->logdie("input file $options{'input'} does not exist")
            unless(-e "$options{input}.gz");
        $options{'input'}.='.gz';
    }
    
    ## make user an output file was passed
    if (! $options{'output'}) { $logger->logdie("output option required!") }

    if($options{'query_file_path'}) {
        open(IN, "< $options{'query_file_path'}") or 
            $logger->logdie("cannot open $options{query_file_path}");
        while(<IN>) {
            chomp;
            if(/^>(.*)/) {
                $qfDefline = $1;
                last;
            }
        }
        close(IN);
    }

    return 1;
}

sub min {
    my ($num1, $num2) = @_;
    
    if ($num1 < $num2) {
        return $num1;
    } else {
        return $num2;
    }
}

##Adds BSML tags for the case where
##the query sequence returned no hits
sub createAndAddNullResult {
    my %args = @_;
    my $doc = $args{'doc'};
     
    if ( !( $doc->returnBsmlSequenceByIDR( "$args{'query_name'}")) ){
        my $seq = $doc->createAndAddSequence( "$args{'query_name'}", "$args{'query_name'}", '', 'na', $args{'class'} );
        $doc->createAndAddSeqDataImport($seq, 'fasta', $options{'query_file_path'}, '', $args{'query_name'});
        $seq->addBsmlLink('analysis', '#aat_aa_analysis', 'input_of');
    }
}
