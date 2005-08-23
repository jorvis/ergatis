#!/usr/local/bin/perl

=head1  NAME 

repeatmasker2bsml.pl - convert RepeatMasker output to BSML

=head1 SYNOPSIS

USAGE: repeatmasker2bsml.pl --input=/path/to/somefile.out --output=/path/to/output.bsml

=head1 OPTIONS

B<--input,-i> 
    Input .out file from a RepeatMasker search.

B<--debug,-d> 
    Debug level.  Use a large number to turn on verbose debugging. 

B<--log,-l> 
    Log file

B<--output,-o> 
    Output BSML file (will be created, must not exist)

B<--help,-h> 
    This help message

=head1   DESCRIPTION

This script is used to convert the output from a RepeatMasker search into BSML.

=head1 INPUT

RepeatMasker can be run using multiple input sequences simultaneously, and this
script supports parsing single or multiple input result sets.  It generates a few
output files, but the tab-delimited ".out" file is used here.  A usual RepeatMasker
.out file looks like (wide-window):

       SW  perc perc perc  query     position in query               matching       repeat         position in  repeat
    score  div. del. ins.  sequence     begin      end     (left)   repeat         class/family   begin  end (left)   ID

      254  22.7 10.5  0.8  id51595      30597    30729 (19616362) +  (AGCTG)n       Simple_repeat      5  150    (0)    1  
      204  25.8  0.0  6.6  id51595      31229    31365 (19615726) +  (CGGA)n        Simple_repeat      1  128    (0)    2  
      237  23.5 10.5  0.8  id51595      33531    33663 (19613428) +  (AGCTG)n       Simple_repeat      5  150    (0)    3  
      204  25.8  0.0  6.6  id51595      34163    34299 (19612792) +  (CGGA)n        Simple_repeat      1  128    (0)    4  
      252   0.0  0.0  0.0  id51595      37646    37673 (19609418) +  (A)n           Simple_repeat      1   28    (0)    5  
       21   0.0  0.0  0.0  id51595      49671    49691 (19597400) +  AT_rich        Low_complexity     1   21    (0)    6  
       21   0.0  0.0  0.0  id51595      52652    52672 (19594419) +  AT_rich        Low_complexity     1   21    (0)    7  

The header rows should be ignored by this script.

You define the input file using the --input option.  This file does not need any
special file extension.

=head1 OUTPUT

After parsing the input file, a file specified by the --output option is created.  This script
will fail if it already exists.  The file is created, and temporary IDs are created for
each result element.  They are only unique to the document, and will need to be replaced
before any database insertion.

Base positions from the input file are renumbered so that positions start at zero.  The
current output elements from RepeatFinder that are not represented in the BSML file are:

    SW score
    perc_div
    perc_del
    perc_ins
    matching_repeat ?
    repeat_class/family

These need to be included later.

=head1 CONTACT

    Joshua Orvis
    jorvis@tigr.org

=cut

use strict;
use Log::Log4perl qw(get_logger);
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use BSML::BsmlBuilder;
use BSML::BsmlReader;
use BSML::BsmlParserTwig;
use BSML::BsmlRepository;
use Pod::Usage;
use Workflow::Logger;

my %options = ();
my $results = GetOptions (\%options, 
			  'input|i=s',
              'output|o=s',
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

## open the input file for parsing
open (my $ifh, $options{'input'}) || $logger->logdie("can't open input file for reading");

my %data;
while (<$ifh>) {
    my @cols = split;
    
    #check whitespace, no warn
    next if ( /^\s*$/ );
    
    ## make sure we don't parse the repeatmasker output header lines
    next if ( /^sw.*position.*repeat/i ||
              /^score.*class.*id/i);

    ## there should be 15 elements in cols, unless we have an unrecognized format.
    unless (scalar @cols == 15) {
        $logger->error("the following RepeatMasker line was not recognized and could not be parsed:\n$_\n") if ($logger->is_error);
        next;
    }
    
    ## add this data row to this sequence
    push( @{$data{$cols[4]}}, \@cols );
}

## loop through each of the matches that we found
for my $seqid (keys %data) {
    my $seq = $doc->createAndAddSequence($seqid);
       $seq->addBsmlLink('analysis', 'repeatmasker_analysis');
    my $ft  = $doc->createAndAddFeatureTable($seq);
    my $fg;
    
    ## loop through each array reference of this key
    my $gene;
    my $repeat;
    my @elements;
    foreach my $arr ( @{$data{$seqid}} ) {
        ## add the repeat
        $repeat = $doc->createAndAddFeature($ft, &fake_id('rep'), '', 'repeat_region');
        $repeat->addBsmlLink('analysis', 'repeatmasker_analysis');
        
        ## add the location of the repeat (all given by RepeatMasker as coords on the forward strand)
        ## 1 is subtracted from each position to give interbase numbering
        $repeat->addBsmlIntervalLoc( --$$arr[5], --$$arr[6], 0);
        
        ## SO terms for these repeats need to be added??
    }
}

## add the analysis element
$doc->createAndAddAnalysis(
                            id => 'repeatmasker_analysis',
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
    
    return 1;
}

sub fake_id {
    my $so_type = shift;

    ## this will be used by the id replacement software to find ids to replace.
    my $temp_prefix = 'ir';
    
    return "$temp_prefix.$so_type." . $next_id++;
}
