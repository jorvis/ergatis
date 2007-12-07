#!/usr/bin/perl
use warnings;
use strict;
$|++;

BEGIN{foreach (@INC) {s/\/usr\/local\/packages/\/local\/platform/}};
use lib (@INC,$ENV{"PERL_MOD_DIR"});
no lib "$ENV{PERL_MOD_DIR}/i686-linux";
no lib ".";

=head1 NAME

seg2bsml.pl - convert seg output to bsml.

=head1 SYNOPSIS

    USAGE: seg2bsml.pl
                --input|-i seg.raw
                --seg_input_file|-s model.fsa
                --output|-o seg.bsml
                --project|-p foo
                --id_repository|-r /path/to/repository
              [ --log|-l seg2bsml.log
                --help
                --debug  (currently unimplemented)
              ]


=head1 OPTIONS

B<--input,-i>
    - The raw seg output (from a seg run using the -l option)

B<--seg_input_file, -s>
    - The file used as input for seg

B<--help,-h>
    - Display this help

=head1  DESCRIPTION

    This program converts raw seg output into bsml.  Currently, the encoding is of
    type repeat_region, but when a satisfactory term is added to the ontology this
    script should be updated to more accurately represent the masked regions.

=head1  INPUT

    It is expected that the input be in the format that seg uses for its output when
    invoked with the -l option.  An example of the format expected:


>PRIO_HUMAN(50-94) complexity=1.92 (12/2.20/2.50)
ppqggggwgqphgggwgqphgggwgqphgggwgqphgggwgqggg

>PRIO_HUMAN(113-135) complexity=2.47 (12/2.20/2.50)
agaaaagavvgglggymlgsams

>PRIO_HUMAN(188-201) complexity=2.26 (12/2.20/2.50)
tvttttkgenftet


=head1  OUTPUT

    The output is a bsml file with the masked regions described as 'repeat-regions'.

=head1  CONTACT

    Jason Inman
    jinman@tigr.org

=begin comment
    ## legal values for status are active, inactive, hidden, unstable
    status: active
    keywords: keywords to search for your script here
=end comment

=cut

use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use Pod::Usage;
use Ergatis::IdGenerator;
use Ergatis::Logger;
use BSML::BsmlRepository;
use BSML::BsmlBuilder;
use BSML::BsmlParserTwig;

my %opts = ();
my $input_file;
my $output_file;
my $seg_input_file;
my $id_repository;
my $project;

my $results = GetOptions(\%opts,
                        'input|i=s',
                        'seg_input_file|s=s',
                        'output|o=s',
                        'project|p=s',
                        'log|l=s',
                        'id_repository|r=s',
                        'command_id=s',     # passed by workflow
                        'logconf=s',        # passed by workflow (not used)
                        'debug=s',
                        'help|h',
                        ) || pod2usage();

my $logfile = $opts{'log'} || Ergatis::Logger::get_default_logfilename();
my $logger = new Ergatis::Logger('LOG_FILE'=>$logfile, 'LOG_LEVEL'=>$opts{'debug'});
$logger = $logger->get_logger();

# display documentation
if( $opts{'help'} ){
    pod2usage( {-exitval=>0, -verbose => 2, -output => \*STDOUT} );
}

## make sure all passed options are peachy
&check_parameters(\%opts);

## we want a new doc
my $doc = new BSML::BsmlBuilder();

## we're going to generate ids
my $idcreator = new Ergatis::IdGenerator( id_repository => $id_repository );

## open the input file for parsing
open (my $ifh, $input_file) || $logger->logdie("can't open input file for reading");

## parse the data from the input file
my @data = ();
my $seqid;

while (<$ifh>) {

    next unless (/^>(\S+)\((\d+)-(\d+)\) complexity=(\d+\.\d+ \([\d\.\/]+\))\n/);
    my %row;

    if ($seqid) {

        ## Normally, seg is fine for this, but for our bsml, we don't like it:
        $logger->logdie("Multiple input seqs detected!") if ($seqid ne $1);

    } else {
        $seqid = $1;
    }

    $row{'start'}   = $2;
    $row{'end'}     = $3;
    $row{'complexity'} = $4;

    push @data, \%row;
        
}

## Begin putting seq data into the bsml doc
my $seq = $doc->createAndAddSequence($seqid, undef, '', 'aa', 'polypeptide');
$seq->addBsmlLink('analysis','#seg_analysis','input_of');

## seq data import:
$doc->createAndAddSeqDataImport( $seq, 'fasta', $seg_input_file, '', $seqid );

## feature table:
my $ft = $doc->createAndAddFeatureTable($seq);
my $mask;
foreach my $row (@data) {

    # create and add the mask
    ## Remember to change this away from repeat_region to masked region or whatever the 
    ## ontology stores it as eventually
    my $id = $idcreator->next_id( project => $project, type => 'repeat_region' );
    $mask = $doc->createAndAddFeature($ft,$id, '','repeat_region');
    $mask->addBsmlLink('analysis','#seg_analysis','computed_by');

    ## make it interbase numbering:
    $$row{'start'}--;
    $$row{'end'}--;
    $mask->addBsmlIntervalLoc( $$row{'start'}, $$row{'end'}, 0);

    ## add the properties (complexity) to the masked region.
    $doc->createAndAddBsmlAttributes( $mask, 'complexity', $$row{'complexity'} );


}

## add the analysis:
$doc->createAndAddAnalysis( id => 'seg_analysis', sourcename => $output_file );

## write it out:
$doc->write($output_file);

## Bob's yer uncle, mate!
exit(0);

sub check_parameters {

    my $error = '';

    # input file checking
    if (exists $opts{'input'}) {
        $input_file = $opts{'input'};
        ## make sure input file exists
        $error .= "Input file $input_file does not exist.\n" unless (-e $input_file);        
    } else {
        $error .= "No input file specified.\n";
    } 

    # output file checking
    if (exists $opts{'output'}) {
        $output_file = $opts{'output'};
    } else {
        $error .= "No output file specified.\n";
    }
    
    # housekeeping stuff
    if (exists $opts{'project'}) {
        $project = $opts{'project'};
    } else {
        $error .= "No project specified.\n";
    }
    if (exists $opts{'id_repository'}) {
        $id_repository = $opts{'id_repository'};
        unless (-e "$id_repository/valid_id_repository") {
            $error .= "$id_repository does not appear to be a valid id repository.\n";
        }
    if (exists $opts{'seg_input_file'}) {
        $seg_input_file = $opts{'seg_input_file'};
        unless (-e $seg_input_file) {
            $error .= "Can't find (seg input) $seg_input_file.\n";
        }
    } else {
        $error .= "No seg_input_file specified.\n";
    }

    } else {
        $error .= "No id_repository specified.\n";
    }

    if ($error) {
        $logger->logdie("$error"); 
    }

}

