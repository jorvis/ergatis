#!/usr/bin/perl

eval 'exec /local/packages/perl-5.8.8/bin/perl  -S $0 ${1+"$@"}'
    if 0; # not running under some shell
use lib (@INC,$ENV{"PERL_MOD_DIR"});
no lib "$ENV{PERL_MOD_DIR}/i686-linux";
no lib ".";

=head1  NAME 

mreps2bsml.pl - convert mreps output to BSML

=head1 SYNOPSIS

USAGE: mreps2bsml.pl 
            --input=/path/to/somefile.dat 
            --output=/path/to/output.bsml
          [ --project=aa1 ]

=head1 OPTIONS

B<--input,-i> 
    Input .raw file from an mreps search.

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

This script is used to convert the output from an mreps search into BSML.

=head1 INPUT

mreps output is to stdout.  mreps2bsml works on the captured output of mreps,
an example of which appears below:

 *****************************************************************************  
 *                              mreps 2.5                                    *  
 *                                                                           *  
 *                Finding tandem repeats in DNA sequences                    *  
 *                                                                           *  
 *                      http://www.loria.fr/mreps/                           *  
 *****************************************************************************  

Processing sequence 'ThisIsAHeader'

* Processing window [1 : 2216158] *

Warning: symbols N in your sequence has been replaced randomly
by {A,C,G,T}. This might have created artefact repetitions.

   from   ->       to  :         size    <per.>  [exp.]          err-rate       sequence
 ---------------------------------------------------------------------------------------------
     158  ->       215 :         58      <1>     [58.00]         0.333          T A A A T C T T T A G C C T A A T T T T A T T T T A T T T T G T T T T T A T T A T T T A T T T A A A A C T T T T A A
     171  ->       205 :         35      <4>     [8.75]          0.226          TAAT TTTA TTTT ATTT TGTT TTTA TTAT TTAT TTA
     171  ->       196 :         26      <5>     [5.20]          0.190          TAATT TTATT TTATT TTGTT TTTAT T
     395  ->       419 :         25      <4>     [6.25]          0.238          TGTA TATT TAAT TAAT TAAT TATT A

The header rows should be ignored by this script.

Not all columns described in the header are seperated by tabs.  The columns as described
by the headers:
1) 1-based start
2) 1-based end
3) size of the entire repeat region
4) period - size of the repeat unit
5) exponent - number of repet units in the region
6) error rate - lower means more significant
7) repeat sequence, with spaces seperating units of the repeat

You define the input file using the --input option.  This file does not need any
special file extension.  

=head1 OUTPUT

After parsing the input file, a file specified by the --output option is created.  This script
will fail if it already exists.  The file is created, and temporary IDs are created for
each result element. 

Base positions from the input file are renumbered so that positions start at zero.  

=head1 CONTACT

    Jason Inman
    jinman@tigr.org

=cut

use strict;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use Pod::Usage;
BEGIN {
use Ergatis::Logger;
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

## we want to create ids unique to this document, which will be replaced later.  they must
##  contain the prefix that will be used to look up a real id, such as ath1.gen.15
my $next_id = 1;

## we want a new doc
my $doc = new BSML::BsmlBuilder();

## we're going to generate ids
my $idcreator = new Papyrus::TempIdCreator();

## open the input file for parsing
open (my $ifh, $options{'input'}) || $logger->logdie("can't open input file for reading");

my @data;
my $qry_id;
my $parameters;
while (<$ifh>) {

    ## Remove any preceding whitespace
    chomp;
    $_ =~ s/^\s*//;

    # Capture the header if this is that line:
    if (/Processing sequence '(.*)'/) {
        $qry_id = $1;
    }

    ## only the data lines in the output file start with numbers
    if ( /^\d/ ) {
        # Due to the tricky formating of the output, let's force the split 
        # to create 9 elements in the array.  We'll just ignore certain columns.
        my @cols = split(/\s+/,$_,9);

        ## add this data row to this sequence
        push( @data, \@cols );

    }

}

## set up the bsml for the sequence
my $seq = $doc->createAndAddSequence($qry_id, undef, '', 'dna', 'assembly');
   $seq->addBsmlLink('analysis', '#mreps_analysis', 'input_of');
my $ft  = $doc->createAndAddFeatureTable($seq);
my $fg;

## loop through each array reference now
my $repeat;
my @elements;
foreach my $arr ( @data ) {
    ## grab an ID
    my $new_id = $idcreator->new_id( db     => $options{project},
                                                      so_type => 'tandem_repeat',
                                                      prefix => $options{command_id}
                                                    );
    
    ## add the repeat
    $repeat = $doc->createAndAddFeature($ft, $new_id, '', $idcreator->so_used('tandem_repeat') );
    $repeat->addBsmlLink('analysis', '#mreps_analysis', 'computed_by');
        
    ## add the location of the repeat (all given by mreps as coords are on the forward strand)
    ## 1 is subtracted from each position to give interbase numbering
    $repeat->addBsmlIntervalLoc( --$$arr[0], --$$arr[2], 0);
   
    ## Adjust the vlaues in period and exponent cells to take them out of the silly
    ## brackets.  Also strip out the spaces since we know the periodicity and exponent.
    $$arr[5] =~ s/<(.*)>/$1/;
    $$arr[6] =~ s/\[(.*)\]/$1/;
    $$arr[8] =~ tr/ //d;
    
 
    ## SO terms for these repeats need to be added as Attributes
    $doc->createAndAddBsmlAttributes( $repeat, 'period_size',       $$arr[5],
                                               'exponent',          $$arr[6],
                                               'consensus_size',    $$arr[4],
                                               'err-rate',          $$arr[7],
                                               'repeat_sequence',   $$arr[8],
                                    );
}


## add the analysis element
my $analysis = $doc->createAndAddAnalysis(
                            id => 'mreps_analysis',
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

