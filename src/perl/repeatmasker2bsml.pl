#!/usr/local/bin/perl
BEGIN{foreach (@INC) {s/\/usr\/local\/packages/\/local\/platform/}};
use lib (@INC,$ENV{"PERL_MOD_DIR"});
no lib "$ENV{PERL_MOD_DIR}/i686-linux";
no lib ".";

=head1  NAME 

repeatmasker2bsml.pl - convert RepeatMasker output to BSML

=head1 SYNOPSIS

USAGE: repeatmasker2bsml.pl 
        --input=/path/to/somefile.out 
        --output=/path/to/output.bsml
      [ --project=aa1 ]

=head1 OPTIONS

B<--input,-i> 
    Input .out file from a RepeatMasker search.

B<--debug,-d> 
    Debug level.  Use a large number to turn on verbose debugging. 

B<--log,-l> 
    Log file

B<--output,-o> 
    Output BSML file (will be created, or overwritten)

B<--project,-p> 
    Project ID.  Used in creating feature ids.  Defaults to 'unknown' if
    not passed.

B<--id_repository,-r> 
    Required for creating feature identifiers.  Each project should have
    its own id_repository directory - use the full path to it here.

B<--help,-h> 
    This help message

=head1   DESCRIPTION

This script is used to convert the output from a RepeatMasker search into BSML.

=head1 INPUT

RepeatMasker can be run using multiple input sequences simultaneously, and this
script supports parsing single or multiple input result sets.  It generates a few
output files, but the space-delimited ".out" file is used here.  A usual RepeatMasker
.out file looks like (wide-window):

       SW  perc perc perc  query                position in query           matching       repeat             position in repeat
    score  div. del. ins.  sequence               begin      end   (left)   repeat         class/family           begin  end   (left)

       219  6.9  0.0  0.0 sma1.assembly.30997     11239    11267  (23517) +  (TA)n          Simple_repeat         1      29     (151)
       572 24.2  4.0  5.7 sma1.assembly.30997     11493    11720  (23045) C  R=109                             (226)    1736    1513 *
      1124 17.9  0.0  2.2 sma1.assembly.30997     11539    11766  (22999) C  R=86                                (1)     511     289 *

The * at the end of some lines indicates an overlap with other predictions.

You define the input file using the --input option.  This file does not need any
special file extension.

=head1 OUTPUT

After parsing the input file, a file specified by the --output option is created.  
The file is created, and temporary IDs are created for each result element.  They are 
globally unique, but will need to be replaced before any database insertion.

Base positions from the input file are renumbered so that positions start at zero.  The
current output elements from RepeatFinder that are not represented in the BSML file are:

    perc_div
    perc_del
    perc_ins

We can add these later if anyone thinks they are necessary.

=head1 CONTACT

    Joshua Orvis
    jorvis@tigr.org

=cut

use strict;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use Pod::Usage;
use Workflow::Logger;
use Workflow::IdGenerator;
use BSML::BsmlRepository;
use BSML::BsmlBuilder;
use BSML::BsmlParserTwig;

my %options = ();
my $results = GetOptions (\%options, 
              'input|i=s',
              'output|o=s',
              'project|p=s',
              'log|l=s',
              'id_repository|r=s',
              'command_id=s',       ## passed by workflow
              'logconf=s',          ## passed by workflow (not used)
              'debug=s',
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

## we want a new doc
my $doc = new BSML::BsmlBuilder();

## we're going to generate ids
my $idcreator = new Workflow::IdGenerator( id_repository => $options{id_repository} );

## open the input file for parsing
open (my $ifh, $options{'input'}) || $logger->logdie("can't open input file for reading");

my %data;

my $linectr=0;
my $datalinectr=0;

my $old_version=0;

while (<$ifh>) {


    $linectr++;

    chomp;
    #check whitespace, no warn
    next if ( /^\s*$/ );
    
    ## make sure we don't parse the repeatmasker output header lines
    ## (I think these don't exist anymore in current versions)
    if (( /^\s*sw.*position.*repeat/i) || (/^\s*score.*class.*id/i)) {

        ## set flag which indicates script is parsing older version of repeatmasker output
        $old_version = 1;

        next;
    }

    ## let's count the data lines
    $datalinectr++;

    ## skip lines that warn about lack of matches:
    next if ( /were no repetitive sequences detected/i );

    ## if the line ends with an asterisk, remove it.  it only causes
    ## column ambiguity.
    s/(.+)\s+\*\s*$/$1/;

    my @cols = split;

    ## get the number of columns parsed
    my $colctr = scalar(@cols);

    
    if (($old_version) && ($colctr == 15)){

        ## if we are parsing the old version of input file and the number of columns is fifteen with last column being an integer (ID value)- remove that last column.  it only causes
        ## column ambiguity.
        s/(.+)\s+\*\s*$/$1/;

        $colctr--;
    }



    ## if there are 13 columns, the repeat family must have been missing and we need to adjust.
    if ( $colctr == 13 ) {
    ( $cols[10], $cols[11], $cols[12], $cols[13] ) = ( '', $cols[10], $cols[11], $cols[12] );
    
        push( @{$data{$cols[4]}}, \@cols );    
    
    ## if there are 14 columns, just add.
    } elsif ( $colctr ==14 ) {
        push( @{$data{$cols[4]}}, \@cols );
    
    ## else we have an unrecognized row.
    } else {
        $logger->logdie("Could not parse file '$options{'input'}' at file line '$linectr' data line '$datalinectr'.  Number of columns '$colctr'.  The following RepeatMasker line was not recognized and could not be parsed:\n$_\n");
    }


}

## loop through each of the matches that we found
for my $seqid (keys %data) {
    my $seq = $doc->createAndAddSequence($seqid, undef, '', 'dna', 'assembly');
       $seq->addBsmlLink('analysis', '#repeatmasker_analysis', 'input_of');
    my $ft  = $doc->createAndAddFeatureTable($seq);
    my $fg;
    
    ## loop through each array reference of this key
    my $gene;
    my $repeat;
    my @elements;
    foreach my $arr ( @{$data{$seqid}} ) {
        ## add the repeat
        my $id = $idcreator->next_id( project => $options{project}, type => 'repeat_region' );
        $repeat = $doc->createAndAddFeature($ft, $id, '', 'repeat_region');
        $repeat->addBsmlLink('analysis', '#repeatmasker_analysis', 'computed_by');
        
        ## add the location of the repeat
        ## 1 is subtracted from each position to give interbase numbering
        if ($$arr[8] eq '+') {
            $repeat->addBsmlIntervalLoc( --$$arr[5], $$arr[6], 0);
        } elsif ($$arr[8] eq 'C') {
            $repeat->addBsmlIntervalLoc( --$$arr[5], $$arr[6], 1);
        } else {
            $logger->logdie("expected '+' or 'C' in column 9, but found: $$arr[8]");
        }
        
        ## add the properties of the repeat:
        $doc->createAndAddBsmlAttributes( $repeat, 'matching_repeat',   $$arr[9],
                                                   'sw_score',          $$arr[0],
                                        );

        ## add the repeat class, if we have one
        if ($$arr[10]) {
            $doc->createAndAddBsmlAttributes( $repeat, 'repeat_class', $$arr[10] );
        }

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
    
    return 1;
}
