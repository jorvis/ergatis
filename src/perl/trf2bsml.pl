#!/usr/local/bin/perl
use lib (@INC,$ENV{"PERL_MOD_DIR"});
no lib "$ENV{PERL_MOD_DIR}/i686-linux";
no lib ".";

=head1  NAME 

trf2bsml.pl - convert Tandem Repeat Finder (trf) output to BSML

=head1 SYNOPSIS

USAGE: trf2bsml.pl 
            --input=/path/to/somefile.dat 
            --output=/path/to/output.bsml
          [ --project=aa1 ]

=head1 OPTIONS

B<--input,-i> 
    Input .dat file from a RepeatMasker search.

B<--debug,-d> 
    Debug level.  Use a large number to turn on verbose debugging. 

B<--log,-l> 
    Log file

B<--output,-o> 
    Output BSML file (will be created, must not exist)

B<--project,-p>
    Project ID.  Used in creating feature ids.  Defaults to 'unknown' if
    not passed.

B<--help,-h> 
    This help message

=head1   DESCRIPTION

This script is used to convert the output from a RepeatMasker search into BSML.

=head1 INPUT

trf can be run using a query file with one or many input sequences, and this
script supports parsing single or multiple input result sets.  trf generates a few
output files, but the space-delimited ".dat" file is used here.  A usual trf
.dat file looks like (wide-window):

    Sequence: cpa1.assem.1

    Parameters: 2 7 7 80 10 50 500

    34540 34614 3 25.0 3 82 15 93 68 1 30 0 0.98 AGA AGAAGAAGAAGAAGAAGAAGAAAACGAAGAAGAAGAAGAAAGAAAAGAAGAAGAAGAGAGAGAAGAAGAAAAAGA
    35458 35499 3 13.7 3 95 5 75 64 0 2 33 1.07 ATA ATAATGAATAATAATAATAATAATAATAATAATAATAATAAT
    37239 37288 3 16.3 3 95 4 91 32 2 0 66 1.03 TTA TTATTATTATTATTATTATTATTATTATTATTATTATTATTATTCATTAT
    48180 48226 3 15.3 3 95 4 85 65 0 2 31 1.04 ATA ATAATGAATAATAATAATAATAATAATAATAATAATAATAATAATAA
    51008 51037 3 9.7 3 92 7 51 33 0 0 66 0.92 ATT ATTATTATTATTATTATTATTATTTATTAT 1   21    (0)    7  

The header rows should be ignored by this script.

You define the input file using the --input option.  This file does not need any
special file extension.

=head1 OUTPUT

After parsing the input file, a file specified by the --output option is created.  This script
will fail if it already exists.  The file is created, and temporary IDs are created for
each result element. 

Base positions from the input file are renumbered so that positions start at zero.  

=head1 CONTACT

    Joshua Orvis
    jorvis@tigr.org

=cut

use strict;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use Pod::Usage;
BEGIN {
use Workflow::Logger;
use BSML::BsmlRepository;
use Papyrus::TempIdCreator;
use BSML::BsmlBuilder;
use BSML::BsmlParserTwig;
}

my %options = ();
my $results = GetOptions (\%options, 
			  'input|i=s',
              'output|o=s',
              'debug|d=s',
              'command_id=s',       ## passed by workflow
              'logconf=s',          ## passed by workflow (not used)
              'project|p=s',
              'log|l=s',
			  'help|h') || pod2usage();

my $logfile = $options{'log'} || Workflow::Logger::get_default_logfilename();
my $logger = new Workflow::Logger('LOG_FILE'=>$logfile,
				  'LOG_LEVEL'=>$options{'debug'});
$logger = $logger->get_logger();

# display documentation
if( $options{'help'} ){
    pod2usage( {-exitval=>0, -verbose => 2, -output => \*STDOUT} );
}

## make sure all passed options are peachy
&check_parameters(\%options);

## we want to create ids unique to this document, which will be replaced later.  they must
##  contain the prefix that will be used to look up a real id, such as ath1.gen.15
my $next_id = 1;

## we want a new doc
my $doc = new BSML::BsmlBuilder();

## we're going to generate ids
my $idcreator = new Papyrus::TempIdCreator();

## open the input file for parsing
open (my $ifh, $options{'input'}) || $logger->logdie("can't open input file for reading");

my %data;
my $qry_id;
my $parameters;
while (<$ifh>) {
    #check whitespace, no warn
    next if ( /^\s*$/ );

    ## only the data lines in the output file start with numbers
    if ( /^\d/ ) {
        my @cols = split;

        ## there should be 15 elements in cols, unless we have an unrecognized format.
        unless (scalar @cols == 15) {
            $logger->error("the following RepeatMasker line was not recognized and could not be parsed:\n$_\n") if ($logger->is_error);
            next;
        }

        ## add this data row to this sequence
        push( @{$data{$qry_id}}, \@cols );

    } elsif ( /^Sequence\:\s*(.+)\s*$/ ) {
        $qry_id = $1;
    
    } elsif ( /^Parameters\:\s*(.+)\s*$/ ) {
        $parameters = $1;
    }
}

## loop through each of the matches that we found
for my $seqid (keys %data) {
    my $seq = $doc->createAndAddSequence($seqid, undef, '', 'dna', 'assembly');
       $seq->addBsmlLink('analysis', '#trf_analysis', 'input_of');
    my $ft  = $doc->createAndAddFeatureTable($seq);
    my $fg;
    
    ## loop through each array reference of this key
    my $repeat;
    my @elements;
    foreach my $arr ( @{$data{$seqid}} ) {
        ## grab an ID
        my $new_id = $idcreator->new_id( db     => $options{project},
                                                          so_type => 'tandem_repeat',
                                                          prefix => $options{command_id}
                                                        );
        
        ## add the repeat
        $repeat = $doc->createAndAddFeature($ft, $new_id, '', $idcreator->so_used('tandem_repeat') );
        $repeat->addBsmlLink('analysis', '#trf_analysis', 'computed_by');
        
        ## add the location of the repeat (all given by RepeatMasker as coords on the forward strand)
        ## 1 is subtracted from each position to give interbase numbering
        $repeat->addBsmlIntervalLoc( --$$arr[0], --$$arr[1], 0);
        
        ## SO terms for these repeats need to be added as Attributes
        $doc->createAndAddBsmlAttributes( $repeat, 'period_size',        $$arr[2],
                                                   'copies_aligned',     $$arr[3],
                                                   'consensus_size',     $$arr[4],
                                                   'percent_identity',   $$arr[5],
                                                   'percent_indels',     $$arr[6],
                                                   'raw_score',          $$arr[7],
                                                   'percent_a',          $$arr[8],
                                                   'percent_c',          $$arr[9],
                                                   'percent_g',          $$arr[10],
                                                   'percent_t',          $$arr[11],
                                                   'entropy_score',      $$arr[12],
                                                   'consensus_text',     $$arr[13],
                                                   'matched_text',       $$arr[14]
                                        );
    }
}

## add the analysis element
my $analysis = $doc->createAndAddAnalysis(
                            id => 'trf_analysis',
                            sourcename => $options{'output'},
                          );

## now write the doc
$doc->write($options{'output'});

exit;

sub check_parameters {
    
    ## make sure input file exists
    if (! -e $options{'input'}) { $logger->logdie("input file $options{'input'} does not exist") }
    
    ## make sure output file doesn't exist yet
    if (-e $options{'output'}) { $logger->logdie("can't create $options{'output'} because it already exists") }
    
    $options{'project'}    = 'unknown' unless ($options{'project'});
    $options{'command_id'} = '0' unless ($options{'command_id'});
    
    return 1;
}

