#!/usr/local/bin/perl
BEGIN{foreach (@INC) {s/\/usr\/local\/packages/\/local\/platform/}};
use lib (@INC,$ENV{"PERL_MOD_DIR"});
no lib "$ENV{PERL_MOD_DIR}/i686-linux";
no lib ".";

=head1  NAME 

glimmer32bsml.pl - convert glimmer3 output to BSML

=head1 SYNOPSIS

USAGE: glimmer32bsml.pl 
        --input=/path/to/glimmer3.output.file 
        --output=/path/to/output.bsml
        --id_respository=/id_repository
       [ --project=aa1 
         --fasta_file=/path/to/fasta/for/Seq-data-import
       ]

=head1 OPTIONS

B<--input,-i> 
    .predict file from a glimmer3 scan.

B<--debug,-d> 
    Debug level.  Use a large number to turn on verbose debugging. 

B<--fasta_file,-f>
    If passed, will create a Seq-data-import element referencing this
    path.

B<--project,-p> 
    Project ID.  Used in creating feature ids.  Defaults to 'unknown' if
    not passed.

B<--log,-l> 
    Log file

B<--output,-o> 
    Output BSML file (will be created, must not exist)

B<--help,-h> 
    This help message

=head1   DESCRIPTION

This script is used to convert the output from a glimmer3 search into BSML.

=head1 INPUT

You define the input file using the --input option.  This file does not need any
special file extension.  The regular output of glimmer3 looks like this:


>zma1.assembly.5808
orf00001   319296      174  +1     3.95
orf00002      430      969  +1     3.91
orf00003     1197     1051  -1    12.79
orf00004     1258     1368  +1     6.95
orf00005     1541     1681  +2     8.58
orf00007     1889     1755  -3    11.14
orf00008     2014     1874  -2     0.37
orf00009     2153     2359  +2    48.95

Where the columns from left to right contain:
orf_id   start_position   end_position   reading_frame   "raw"_score

A few things to note:
- glimmer3 will find genes spanning the end/beginning of the given sequence, as
if dealing with a circular molecule (unless given the -l or --linear option).  These
are dealt with in bsml by converting the start_pos into a coordinate to the left of
the origin (a negative number).
- The raw_score in the .predict file may be different than that found in the 
.details file, due to adjustments made for the PWM and start codon frequency.
This value from the .predict file can therefore be used to make direct comp-
arisons between predictions.
- The ranges given by glimmer3 differ from those given by glimmer2 because they now
include the stop codon.


=head1 OUTPUT

After parsing the input file, a file specified by the --output option is created.  This script
will fail if it already exists.  The file is created, and temporary IDs are created for
each result element.  They are only unique to the document, and will need to be replaced
before any database insertion.

Base positions from the input file are renumbered so that positions start at zero.

=head1 CONTACT

    Jason Inman
    jinman@tigr.org

=cut

use strict;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use Pod::Usage;
BEGIN {
use Workflow::Logger;
use BSML::BsmlRepository;
use Workflow::IdGenerator;
use BSML::BsmlBuilder;
use BSML::BsmlParserTwig;
}

my %options = ();
my $results = GetOptions (\%options, 
              'input|i=s',
              'output|o=s',
              'fasta_file|f=s',
              'id_respository|r=s',
              'project|p=s',
              'log|l=s',
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

## we want to creating ids unique to this document, which will be replaced later.  they must
##  contain the prefix that will be used to look up a real id, such as ath1.exon.15
my $next_id = 1;

## we want a new doc
my $doc = new BSML::BsmlBuilder();

## we're going to generate ids
my $idcreator = new Workflow::IdGenerator( 'id_respository' => $options{'id_repository'} );
#Set a pool size
$idcreator->set_pool_size( 'gene' => 30, 'exon' => 30, 'CDS' => 30 );

## open the input file for parsing
open (my $ifh, $options{'input'}) || $logger->logdie("can't open input file for reading");


## go forward until we get to the sequence name (id)
my $seq_id;
while (<$ifh>) {
    if ( />(.+)\s*.*$/ ) {
        $seq_id = $1;
        last;
    }
}

## make sure we found an id
unless (defined $seq_id) {
    $logger->logdie("unable to pull seq id from $options{'input'}");
}

## create this sequence, an analysis link, and a feature table
my $seq = $doc->createAndAddSequence($seq_id, undef, '', 'dna', 'assembly');
$seq->addBsmlLink('analysis', '#glimmer3_analysis', 'input_of');
my $ft = $doc->createAndAddFeatureTable($seq);

##  also add a link to the fasta file (Seq-data-import) if requested
if ($options{'fasta_file'}) {
    $doc->createAndAddSeqDataImport( $seq, 'fasta', $options{'fasta_file'}, '', $seq_id);
}

## now look for datalines (gene defs)
my (@cols, $fg, $gene);
my ($gene_start, $gene_stop, $gene_dir);

while (<$ifh>) {

    ## matches lines like "orf00002      430      969  +1     3.91"
    if (/^\s*.+\s+\d+\s+\d+\s+.+\s+.+$/) {
        @cols = split();

        ## set the gene start, stop, and direction
        $gene_start = $cols[1] - 1;
        $gene_stop = $cols[2] - 1;
        ($gene_dir   = $cols[3]) =~ s/(.)\d+/$1/;

        # Adjust the start and stop appropriately for strandedness and whether
        # the origin is being spanned.
        if ($gene_dir eq '-') {
            my $tmp = $gene_start;
            $gene_start = $gene_stop;
            $gene_stop = $tmp;
        }

        # If the gene spans the origin, adjust accordingly...
        if ($gene_stop < $gene_start) {

            ## Get the length.  (Determine the name of the .details file,
            # grep it for "Sequence Length = " and then parse the number
            # of bases from that result.)
            (my $details_file = $options{'input'}) =~ s/predict/detail/;
            die "Can't find details file!\n" unless (stat $details_file);
            my $grep_cmd = "grep 'Sequence length' $details_file";
            my $length = `$grep_cmd`;
            $length =~ s/Sequence length = //;

            # subtract the length from the start to get the coordinate in
            # the opposite direction from the origin.
            $gene_start-= $length;

        }

        $gene = $doc->createAndAddFeature($ft, 
                                          $idcreator->next_id( 'type' => 'gene',
                                                               'project' => $options{'project'}),
                                          '', 'gene'
                                         );
        my $cds = $doc->createAndAddFeature($ft,
                                            $idcreator->next_id( 'type' => 'CDS',
                                                                 'project' => $options{'project'}),
                                         '', 'CDS'
                                         );
        my $exon = $doc->createAndAddFeature($ft,
                                             $idcreator->next_id( 'type' => 'exon',
                                                                  'project' => $options{'project'} ),
                                             '', 'exon'
                                             );
        $fg = $doc->createAndAddFeatureGroup( $seq, '', $gene->returnattr('id') );
        $fg->addBsmlFeatureGroupMember( $gene->returnattr('id'), $gene->returnattr('class') );
        $fg->addBsmlFeatureGroupMember( $cds->returnattr('id'), $cds->returnattr('class') );
        $fg->addBsmlFeatureGroupMember( $exon->returnattr('id'), $exon->returnattr('class') );
        $gene->addBsmlLink('analysis', '#glimmer3_analysis', 'computed_by');
        $cds->addBsmlLink('analysis', '#glimmer3_analysis', 'computed_by');
        $exon->addBsmlLink('analysis', '#glimmer3_analysis', 'computed_by');

        #We may now add the interval loc
        &add_interval_loc($gene, $gene_start, $gene_stop, $gene_dir);
        &add_interval_loc($cds, $gene_start, $gene_stop, $gene_dir);
        &add_interval_loc($exon, $gene_start, $gene_stop, $gene_dir);
        

    }

}

## add the analysis element
$doc->createAndAddAnalysis(
                            id => 'glimmer3_analysis',
                            sourcename => $options{'output'},
                          );

## now write the doc
$doc->write($options{'output'});

exit;

sub add_interval_loc {
    my ($feat, $n, $m, $dir) = @_;
    
    ## was it found on the forward or reverse strand?
    if ($dir eq '+') {
        $feat->addBsmlIntervalLoc($n, $m, 0);
    } elsif ($dir eq '-') {
        $feat->addBsmlIntervalLoc($n, $m, 1);
    } else {
        die "Unknown strand direction found!\n";
    }

}

sub check_parameters {
    
    ## make sure input file exists
    if (! -e $options{'input'}) { $logger->logdie("input file $options{'input'} does not exist") }
    
    ## make sure output file doesn't exist yet
    if (-e $options{'output'}) { $logger->logdie("can't create $options{'output'} because it already exists") }
    
    $options{'fasta_file'} = '' unless ($options{'fasta_file'});
    $options{'project'}    = 'unknown' unless ($options{'project'});
    $options{'command_id'} = '' unless ($options{'command_id'});
    
    return 1;
}

