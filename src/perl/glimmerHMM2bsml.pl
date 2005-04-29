#!/usr/local/bin/perl

=head1  NAME 

glimmerHMM2bsml.pl - convert glimmerHMM output to BSML

=head1 SYNOPSIS

USAGE: glimmerHMM2bsml.pl 
        --input=/path/to/glimmerHMM.output.file 
        --output=/path/to/output.bsml
       [ --project=aa1 ]

=head1 OPTIONS

B<--input,-i> 
    Input file file from a glimmerHMM scan.

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

This script is used to convert the output from a glimmerHMM search into BSML.

=head1 INPUT

You define the input file using the --input option.  This file does not need any
special file extension.  The regular output of glimmerHMM looks like this:

    GlimmerHMM
    Sequence name: aa1.assembly.15595
    Sequence length: 73904 bp
    Sequence C+G%: 39.976185


    Predicted genes/exons

    Gene Exon Strand  Exon            Exon Range      Exon
       #    #         Type                           Length

       1    1  -  Initial        2522       2729      208

       2    1  +  Single         4309       4566      258

       3    1  +  Initial        6567       6594       28
       3    2  +  Internal      10892      10912       21
       3    3  +  Terminal      15171      15343      173

       4    1  +  Initial       15371      17356     1986
       4    2  +  Terminal      23748      23900      153

       5    1  -  Terminal      26294      26316       23
       5    2  -  Initial       32039      32054       16

Each element in BSML requires an id, and the Sequence stub created uses as its ID the
value after the "Sequence name: " label above (up to the first whitespace.)

All other header lines in the file are ignored.

=head1 OUTPUT

After parsing the input file, a file specified by the --output option is created.  This script
will fail if it already exists.  The file is created, and temporary IDs are created for
each result element.  They are only unique to the document, and will need to be replaced
before any database insertion.

Base positions from the input file are renumbered so that positions start at zero.

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
              'fasta_file|f=s',
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

## open the input file for parsing
open (my $ifh, $options{'input'}) || $logger->logdie("can't open input file for reading");


## go forward until we get to the sequence name (id)
my $seq_id;
while (<$ifh>) {
    if ( /^Sequence name: (.+)\s*.*$/ ) {
        $seq_id = $1;
        last;
    }
}

## make sure we found an id
unless (defined $seq_id) {
    $logger->logdie("unable to pull seq id from $options{'input'}");
}

## create this sequence, an analysis link, and a feature table
my $seq = $doc->createAndAddSequence($seq_id);
   $seq->addBsmlLink('analysis', '#glimmerHMM_analysis');
my $ft = $doc->createAndAddFeatureTable($seq);

##  also add a link to the fasta file (Seq-data-import) if requested
if ($options{'fasta_file'}) {
    $doc->createAndAddSeqDataImport( $seq, 'fasta', $options{'fasta_file'}, '', $seq_id);
}

## now look for datalines (exon defs)
my (@cols, $fg, $gene);
my $last_gene_num = 0;
my ($gene_start, $gene_stop, $gene_dir);
while (<$ifh>) {
    ## matches lines like "1    1  -  Initial        2522       2729      208"
    if (/^\s*\d+\s+\d+\s+.\s+.+\s+\d+\s+\d+\s+\d+$/) {
        @cols = split();
        
        ## if this gene number is different than the last one we need to 
        #   add it and start a new feature group
        if ($cols[0] ne $last_gene_num) {
            ## if last_gene_num > 0 we're must not be doing the first one. add 
            ##  the interval loc of the last gene before we move on.
            if ($last_gene_num > 0) {
                if ($gene_dir eq '+') {
                    &add_interval_loc($gene, $gene_start, $gene_stop);
                } elsif ($gene_dir eq '-') {
                    &add_interval_loc($gene, $gene_stop, $gene_start);
                } else {
                    $logger->logdie("unrecognized gene direction ($gene_dir)");
                }
            }
            
            $gene = $doc->createAndAddFeature($ft, &fake_id('gene'), '', 'gene');
            $fg = $doc->createAndAddFeatureGroup( $seq, '', $gene->returnattr('id') );
            $fg->addBsmlFeatureGroupMember( $gene->returnattr('id'), $gene->returnattr('class') );
            $gene->addBsmlLink('analysis', '#glimmerHMM_analysis');
            $last_gene_num = $cols[0];
            
            ## set the new gene start and direction
            $gene_start = $cols[4] - 1;
            $gene_dir   = $cols[2];
        }
        
        ## set the gene stop (gets overridden each time)
        $gene_stop = $cols[5] - 1;
        
        ## add this exon feature
        my $exon = $doc->createAndAddFeature($ft, &fake_id('exon'), '', 'exon');
           $exon->addBsmlLink('analysis', '#glimmerHMM_analysis');
        $fg->addBsmlFeatureGroupMember( $exon->returnattr('id'), $exon->returnattr('class') );
        
        if ($cols[2] eq '+') {
            &add_interval_loc($exon, $cols[4] - 1, $cols[5] - 1);
        } elsif ($cols[2] eq '-') {
            &add_interval_loc($exon, $cols[5] - 1, $cols[4] - 1);
        } else {
            $logger->logdie("unknown value ($cols[2]) in column 3 of input file");
        }
    }
}

## add the last gene interval loc
if ($gene_dir eq '+') {
    &add_interval_loc($gene, $gene_start, $gene_stop);
} elsif ($gene_dir eq '-') {
    &add_interval_loc($gene, $gene_stop, $gene_start);
} else {
    $logger->logdie("unrecognized gene direction ($gene_dir)");
}

## add the analysis element
$doc->createAndAddAnalysis(
                            id => 'glimmerHMM_analysis',
                            sourcename => $options{'output'},
                          );

## now write the doc
$doc->write($options{'output'});

exit;

sub add_interval_loc {
    my ($feat, $n, $m) = @_;
    
    ## was it found on the forward or reverse strand?
    if ($n <= $m) {
        $feat->addBsmlIntervalLoc($n, $m, 0);
    } else {
        $feat->addBsmlIntervalLoc($m, $n, 1);
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

sub fake_id {
    my $so_type = shift;

    ## this will be used by the id replacement software to find ids to replace.
    return "ir.$options{project}.$so_type.$options{'command_id'}" . $next_id++;
}

