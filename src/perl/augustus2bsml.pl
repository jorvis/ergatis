#!/usr/local/bin/perl

=head1  NAME 

augustus2bsml.pl - convert augustus GFF (modified) output to BSML

=head1 SYNOPSIS

USAGE: augustus2bsml.pl 
        --input=/path/to/augustus.output.file.raw 
        --output=/path/to/output.bsml
      [ --project=aa1 
        --fasta_file=/path/to/somefile.fsa 
      ]

=head1 OPTIONS

B<--input,-i> 
    Input file from an augustus scan.

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
    Output BSML file

B<--help,-h> 
    This help message

=head1   DESCRIPTION

This script is used to convert the output from an augustus search into BSML.

=head1 INPUT

You define the input file using the --input option.  This file does not need any
special file extension.  The regular output of augustus looks like this:

    aa1.assembly.24383	AUGUSTUS	initial	18982	19143	.	+	0	transcript_id "g1.1"; gene_id "g1";
    aa1.assembly.24383	AUGUSTUS	intron	19144	19195	.	+	.	transcript_id "g1.1"; gene_id "g1";
    aa1.assembly.24383	AUGUSTUS	internal	19196	19531	.	+	0	transcript_id "g1.1"; gene_id "g1";
    aa1.assembly.24383	AUGUSTUS	intron	19532	19599	.	+	.	transcript_id "g1.1"; gene_id "g1";
    aa1.assembly.24383	AUGUSTUS	terminal	19600	19854	.	+	0	transcript_id "g1.1"; gene_id "g1";
    aa1.assembly.24383	AUGUSTUS	stop_codon	19852	19854	.	+	0	transcript_id "g1.1"; gene_id "g1";
    aa1.assembly.24383	AUGUSTUS	CDS	18982	19143	.	+	0	transcript_id "g1.1"; gene_id "g1";
    aa1.assembly.24383	AUGUSTUS	CDS	19196	19531	.	+	0	transcript_id "g1.1"; gene_id "g1";
    aa1.assembly.24383	AUGUSTUS	CDS	19600	19851	.	+	0	transcript_id "g1.1"; gene_id "g1";

which has the general format of:

    <seqname> <source> <feature> <start> <end> <score> <strand> <frame> transcript_id "X"; gene_id "Y";

This script currently assumes that the input GFF file only contains the output from
the analysis of a single FASTA sequence.  

=head1 OUTPUT

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
use Papyrus::TempIdCreator;

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

## we want a new doc
my $doc = new BSML::BsmlBuilder();

## we're going to generate ids
my $idcreator = new Papyrus::TempIdCreator();

## open the input file for parsing
open (my $ifh, $options{'input'}) || $logger->logdie("can't open input file for reading");

my $seq_id;
my ($seq, $ft, $fg);
my ($last_group_name, $current_group_name, $current_transcript_id);
my ($thing, $id);

## go through the file
while (<$ifh>) {
    ## skip comment lines
    next if (/^\#/);
    chomp;
    my @cols = split(/\t/);
    
    ## has the sequence been defined yet?  it's in the first column
    ##  this should happen on the first row only
    unless ($seq_id) {
        $seq_id = $cols[0];
        $seq_id =~ s/\s//g;
        
        ## create this sequence, an analysis link, and a feature table
        $seq = $doc->createAndAddSequence($seq_id);
        $seq->addBsmlLink('analysis', '#augustus_analysis');
        $ft = $doc->createAndAddFeatureTable($seq);
        
        ##  also add a link to the fasta file (Seq-data-import) if requested
        if ($options{'fasta_file'}) {
            $doc->createAndAddSeqDataImport( $seq, 'fasta', $options{'fasta_file'}, '', $seq_id);
        }
    }

    #$current_group_name = $cols[8] || '';
    if ($cols[8] =~ /transcript_id \".+\"\; gene_id \"(.+)\"/) {
        $current_group_name = $1;
    } else {
        $logger->logdie("column 9 on line $. in unrecognized format");
    }

    ## if column 9 (group) is defined and is different than the last one we need to
    ##  create a new feature group
    if ($current_group_name && $current_group_name ne $last_group_name) {
        
        ## remember this group name
        $last_group_name = $current_group_name;
        
        ## pull a new gene id
        $current_transcript_id = $idcreator->new_id( db      => $options{project},
                                                     so_type => 'gene',
                                                     prefix  => $options{command_id}
                                                   );
        $fg = $doc->createAndAddFeatureGroup( $seq, '', $current_transcript_id );
    }
    
    ## adjust both positions so that we are numbering from zero
    $cols[3]--;
    $cols[4]--;
    
    ## change the + and - symbols in strand column to 0 and 1, respectively
    if ($cols[6] eq '+') {
        $cols[6] = 0;
    } elsif ($cols[6] eq '-') {
        $cols[6] = 1;
    } else {
        $logger->logdie("unknown value ($cols[6]) in strand column.  expected + or -.");
    }
    

    ## handle each of the types we know about
    if ($cols[2] eq 'initial' || $cols[2] eq 'internal' || $cols[2] eq 'terminal' || $cols[2] eq 'single') {
        &add_feature('exon', $cols[3], $cols[4], $cols[6] );
    } elsif ($cols[2] eq 'intron') {
        &add_feature('intron', $cols[3], $cols[4], $cols[6] );
    } elsif ($cols[2] eq 'CDS') {
        &add_feature('CDS', $cols[3], $cols[4], $cols[6] );
    } elsif ($cols[2] eq 'stop_codon') {
        &add_feature('transcription_end_site', $cols[3], $cols[4], $cols[6] );
    } elsif ($cols[2] eq 'stop_codon') {
        &add_feature('transcription_end_site', $cols[3], $cols[4], $cols[6] );
    } else {
        $logger->logdie("unrecognized type: $cols[2]");
    }

}


## add the analysis element
$doc->createAndAddAnalysis(
                            id => 'augustus_analysis',
                            sourcename => $options{'output'},
                          );

## now write the doc
$doc->write($options{'output'});

exit;

sub add_feature {
    my ($type, $start, $stop, $strand) = @_;
    
    $id = $idcreator->new_id( db => $options{project}, so_type => $type, prefix => $options{command_id} );
    $thing = $doc->createAndAddFeature( $ft, $id, '', $idcreator->so_used($type) );
    $thing->addBsmlLink('analysis', '#augustus_analysis');
    $thing->addBsmlIntervalLoc($start, $stop, $strand);

    $fg->addBsmlFeatureGroupMember( $id, $idcreator->so_used($type) );
    
    ## if type is a primary_transcript we need to add a gene too
    if ($type eq 'primary_transcript') {
        $thing = $doc->createAndAddFeature( $ft, $current_transcript_id, '', $idcreator->so_used('gene') );
        $thing->addBsmlLink('analysis', '#augustus_analysis');
        $thing->addBsmlIntervalLoc($start, $stop, $strand);
        $fg->addBsmlFeatureGroupMember( $current_transcript_id, $idcreator->so_used('gene') );
    }
}

sub check_parameters {
    
    ## make sure input file exists
    if (! -e $options{'input'}) { $logger->logdie("input file $options{'input'} does not exist") }
    
    $options{'fasta_file'} = '' unless ($options{'fasta_file'});
    $options{'project'}    = 'unknown' unless ($options{'project'});
    $options{'command_id'} = '0' unless ($options{'command_id'});
    
    return 1;
}

