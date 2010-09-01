#!/usr/bin/perl

=head1  NAME 

islandpath2bsml.pl - convert islandpath output to BSML

=head1 SYNOPSIS

USAGE: islandpath2bsml.pl 
        --input=/path/to/somefile.out 
        --output=/path/to/output.bsml
      [ --project=aa1 ]

=head1 OPTIONS

B<--input,-i> 
    Input .out file from a Islandpath search.

B<--debug,-d> 
    Debug level.  Use a large number to turn on verbose debugging. 

B<--log,-l> 
    Log file

B<--output,-o> 
    Output BSML file (will be created, or overwritten)

B<--project,-p> 
    Project ID.  Used in creating feature ids.  Defaults to whatever the prefix of the
    input file is (i.e. spntigr4.assembly.1.1 would be spntigr4).

B<--id_repository,-r> 
    Required for creating feature identifiers.  Each project should have
    its own id_repository directory - use the full path to it here.

B<--help,-h> 
    This help message

=head1   DESCRIPTION

This script is used to convert the output from an Islandpath run into BSML.

=head1 INPUT

You define the input file using the --input option.  This file does not need any
special file extension.

=head1 OUTPUT

After parsing the input file, a file specified by the --output option is created.  
The file is created, and temporary IDs are created for each result element.  They are 
globally unique, but will need to be replaced before any database insertion.

=head1 CONTACT

    David Riley
    driley@som.umaryland.edu

=cut

use strict;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use Pod::Usage;
use Ergatis::IdGenerator;
use File::Basename;
use Ergatis::Logger;
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

## we're going to generate ids
my $idcreator = new Ergatis::IdGenerator( id_repository => $options{id_repository} );

## Let's pull the sequence ID out of the filename
my $seq_id = basename($options{'input'}, '.out');
my $prefix = $seq_id;
if($seq_id =~ /([^\.]+)\..*/) {
    $prefix = $1;
}
my $project = $options{project} ?  $options{project} : $prefix;

## first we'll add the sequence element
my $seq = $doc->createAndAddSequence($seq_id, undef, '', 'dna', 'assembly');
$seq->addBsmlLink('analysis', '#islandpath_analysis', 'input_of');
my $ft  = $doc->createAndAddFeatureTable($seq);

## open the input file for parsing
open (my $ifh, $options{'input'}) || $logger->logdie("can't open input file for reading");

my %data;

while (<$ifh>) {

    my ($start, $stop) = split(/\s+/, $_);
    if(defined($start) && defined($stop)) {
        ## add the island
        my $id = $idcreator->next_id( project => $project, type => 'pathogenic_island' );
        my $island = $doc->createAndAddFeature($ft, $id, '', 'pathogenic_island');
        $island->addBsmlLink('analysis', '#islandpath_analysis', 'computed_by');
        
        ## add the location of the repeat
        ## 1 is subtracted from each position to give interbase numbering
        if($start < $stop) }
            $island->addBsmlIntervalLoc( $start, $stop, 0);
        }
        else {
            $island->addBsmlIntervalLoc( $stop, $start, 1);
        }
    }
}

## add the analysis element
$doc->createAndAddAnalysis(
                            id => 'islandpath_analysis',
                            sourcename => dirname($options{'output'}),
                            program => 'islandpath',
                            algorithm => 'islandpath',
                            programversion => 'current'
                          );

## now write the doc
$doc->write($options{'output'});

exit;

sub check_parameters {
    
    ## make sure input file exists
    if (! -e $options{'input'}) { $logger->logdie("input file $options{'input'} does not exist") }
    
    return 1;
}
