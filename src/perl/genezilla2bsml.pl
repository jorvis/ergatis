#!/usr/local/bin/perl

=head1  NAME 

genezilla2bsml.pl - convert genezilla GFF output to BSML

=head1 SYNOPSIS

USAGE: genezilla2bsml.pl 
        --input=/path/to/genezilla.output.file.raw 
        --output=/path/to/output.bsml
      [ --project=aa1 
        --fasta_file=/path/to/somefile.fsa 
      ]

=head1 OPTIONS

B<--input,-i> 
    Input file from an genezilla scan.

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

This script is used to convert the output from an genezilla search into BSML.

=head1 INPUT

You define the input file using the --input option.  This file does not need any
special file extension.  The regular output of genezilla looks like this:

    aa1.assembly.20285	genezilla	poly-A-signal	2654	2660	.	-	.	transgrp=1;
    aa1.assembly.20285	genezilla	final-exon	2748	2774	19.9	-	0	transgrp=1;
    aa1.assembly.20285	genezilla	initial-exon	3431	3958	19.05	-	0	transgrp=1;
    aa1.assembly.20285	genezilla	initial-exon	14359	14389	3.318	+	0	transgrp=2;
    aa1.assembly.20285	genezilla	final-exon	14445	14602	266	+	1	transgrp=2;
    aa1.assembly.20285	genezilla	poly-A-signal	15873	15879	.	-	.	transgrp=3;
    aa1.assembly.20285	genezilla	single-exon	15904	16293	38.33	-	0	transgrp=3;

which has the general format of:

    <seqname> <source> <feature> <start> <end> <score> <strand> <frame> transgrp=N;

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
        $logger->debug("adding seq_id $seq_id") if $logger->is_debug;
        $seq = $doc->createAndAddSequence($seq_id);
        $seq->addBsmlLink('analysis', '#genezilla_analysis');
        $ft = $doc->createAndAddFeatureTable($seq);
        
        ##  also add a link to the fasta file (Seq-data-import) if requested
        if ($options{'fasta_file'}) {
            $doc->createAndAddSeqDataImport( $seq, 'fasta', $options{'fasta_file'}, '', $seq_id);
        }
    }

    if ($cols[8] =~ /transgrp=(\d+)\;/) {
        $current_group_name = $1;
    } else {
        $logger->logdie("unrecognized format in column 8: $cols[7]");
    }

    ## if column 8 (group) is defined and is different than the last one we need to
    ##  create a new feature group
    if ($current_group_name && $current_group_name ne $last_group_name) {
        
        ## remember this group name
        $last_group_name = $current_group_name;
        
        ## pull a new gene id
        $current_transcript_id = $idcreator->new_id( db      => $options{project},
                                                     so_type => 'gene',
                                                     prefix  => $options{command_id}
                                                   );
        $logger->debug("adding new feature group with parent $current_transcript_id") if $logger->is_debug;
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
    if ($cols[2] eq 'initial-exon' || $cols[2] eq 'internal-exon' || $cols[2] eq 'final-exon' || $cols[2] eq 'single-exon') {
        &add_feature('exon', $cols[3], $cols[4], $cols[6] );
    } elsif ($cols[2] eq 'poly-A-signal') {
        &add_feature('polyA_signal_sequence', $cols[3], $cols[4], $cols[6] );
    } elsif ($cols[2] eq 'promoter') {
        &add_feature('promoter', $cols[3], $cols[4], $cols[6] );
    } else {
        $logger->logdie("unrecognized type: $cols[2]");
    }

}


## add the analysis element
$doc->createAndAddAnalysis(
                            id => 'genezilla_analysis',
                            sourcename => $options{'output'},
                          );

## now write the doc
$doc->write($options{'output'});

exit;

sub add_feature {
    my ($type, $start, $stop, $strand) = @_;
    
    $id = $idcreator->new_id( db => $options{project}, so_type => $type, prefix => $options{command_id} );
    $thing = $doc->createAndAddFeature( $ft, $id, '', $idcreator->so_used($type) );
    $thing->addBsmlLink('analysis', '#genezilla_analysis');
    $thing->addBsmlIntervalLoc($start, $stop, $strand);

    $fg->addBsmlFeatureGroupMember( $id, $idcreator->so_used($type) );
    
    ## if type is a primary_transcript we need to add a gene too
    if ($type eq 'primary_transcript') {
        $thing = $doc->createAndAddFeature( $ft, $current_transcript_id, '', $idcreator->so_used('gene') );
        $thing->addBsmlLink('analysis', '#genezilla_analysis');
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

