#!/usr/local/bin/perl
BEGIN{foreach (@INC) {s/\/usr\/local\/packages/\/local\/platform/}};
use lib (@INC,$ENV{"PERL_MOD_DIR"});
no lib "$ENV{PERL_MOD_DIR}/i686-linux";
no lib ".";

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

    sma1.assembly.30903     AUGUSTUS        gene    27119   37211   .       +       .       g1
    sma1.assembly.30903     AUGUSTUS        transcript      27119   37211   .       +       .       g1.t1
    sma1.assembly.30903     AUGUSTUS        start_codon     27119   27121   .       +       0       transcript_id "g1.t1"; gene_id "g1";
    sma1.assembly.30903     AUGUSTUS        initial 27119   27245   .       +       0       transcript_id "g1.t1"; gene_id "g1";
    sma1.assembly.30903     AUGUSTUS        intron  27246   27277   .       +       .       transcript_id "g1.t1"; gene_id "g1";
    sma1.assembly.30903     AUGUSTUS        internal        27278   27462   .       +       2       transcript_id "g1.t1"; gene_id "g1";
    sma1.assembly.30903     AUGUSTUS        intron  27463   27508   .       +       .       transcript_id "g1.t1"; gene_id "g1";
    sma1.assembly.30903     AUGUSTUS        internal        27509   27634   .       +       0       transcript_id "g1.t1"; gene_id "g1";
    sma1.assembly.30903     AUGUSTUS        intron  27635   28144   .       +       .       transcript_id "g1.t1"; gene_id "g1";
    sma1.assembly.30903     AUGUSTUS        internal        28145   28294   .       +       0       transcript_id "g1.t1"; gene_id "g1";
    sma1.assembly.30903     AUGUSTUS        intron  28295   30277   .       +       .       transcript_id "g1.t1"; gene_id "g1";
    sma1.assembly.30903     AUGUSTUS        internal        30278   30442   .       +       0       transcript_id "g1.t1"; gene_id "g1";
    sma1.assembly.30903     AUGUSTUS        intron  30443   31648   .       +       .       transcript_id "g1.t1"; gene_id "g1";
    sma1.assembly.30903     AUGUSTUS        internal        31649   31737   .       +       0       transcript_id "g1.t1"; gene_id "g1";
    sma1.assembly.30903     AUGUSTUS        intron  31738   36316   .       +       .       transcript_id "g1.t1"; gene_id "g1";
    sma1.assembly.30903     AUGUSTUS        internal        36317   36669   .       +       1       transcript_id "g1.t1"; gene_id "g1";
    sma1.assembly.30903     AUGUSTUS        intron  36670   36741   .       +       .       transcript_id "g1.t1"; gene_id "g1";
    sma1.assembly.30903     AUGUSTUS        terminal        36742   37211   .       +       2       transcript_id "g1.t1"; gene_id "g1";
    sma1.assembly.30903     AUGUSTUS        stop_codon      37209   37211   .       +       0       transcript_id "g1.t1"; gene_id "g1";
    sma1.assembly.30903     AUGUSTUS        CDS     27119   27245   .       +       0       transcript_id "g1.t1"; gene_id "g1";
    sma1.assembly.30903     AUGUSTUS        CDS     27278   27462   .       +       2       transcript_id "g1.t1"; gene_id "g1";
    sma1.assembly.30903     AUGUSTUS        CDS     27509   27634   .       +       0       transcript_id "g1.t1"; gene_id "g1";
    sma1.assembly.30903     AUGUSTUS        CDS     28145   28294   .       +       0       transcript_id "g1.t1"; gene_id "g1";
    sma1.assembly.30903     AUGUSTUS        CDS     30278   30442   .       +       0       transcript_id "g1.t1"; gene_id "g1";
    sma1.assembly.30903     AUGUSTUS        CDS     31649   31737   .       +       0       transcript_id "g1.t1"; gene_id "g1";
    sma1.assembly.30903     AUGUSTUS        CDS     36317   36669   .       +       1       transcript_id "g1.t1"; gene_id "g1";
    sma1.assembly.30903     AUGUSTUS        CDS     36742   37211   .       +       2       transcript_id "g1.t1"; gene_id "g1";

which has the general format of:

    <seqname> <source> <feature> <start> <end> <score> <strand> <frame> transcript_id "X"; gene_id "Y";

This script currently assumes that the input GFF file only contains the output from
the analysis of a single FASTA sequence.  

=head1 OUTPUT

Base positions from the input file are renumbered so that positions start at zero. Also,
the stop coordinates of each terminal CDS is extended 3bp to include the stop codon.

=head1 CONTACT

    Joshua Orvis
    jorvis@tigr.org

=cut

use strict;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use Pod::Usage;
use Workflow::Logger;
use BSML::BsmlRepository;
use Papyrus::TempIdCreator;
use BSML::BsmlBuilder;
use BSML::BsmlParserTwig;

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

## different versions of augustus have treated stop codons differently.  only some include
##  the stop codon as part of the sequence range.  toggle this to 1 if you want the terminal
##  exon predictions extended 3bp.  this should only be used for those versions of augustus
##  that don't include the stop codon.
my $extend_terminal_exon = 0;

## open the input file for parsing
open (my $ifh, $options{'input'}) || $logger->logdie("can't open input file for reading");

my $seq_id;
my ($seq, $ft, $fg);
my ($last_group_name, $current_group_name, $current_transcript_id);
my ($thing, $id);
## we have to hold each cds in an array so that we can add 3 bases to the last
##  one (augustus doesn't include the stop codon within the cds)
my @cds;

## go through the file
while (<$ifh>) {
    ## skip comment lines
    next if (/^\#/);
    chomp;
    my @cols = split(/\t/);
    
    ## skip the 'gene' and 'transcript' rows, since they are handled by the
    ##  explicit groupings in the last column.
    next if ($cols[2] eq 'gene' || $cols[2] eq 'transcript');
    
    ## has the sequence been defined yet?  it's in the first column
    ##  this should happen on the first row only
    unless ($seq_id) {
        $seq_id = $cols[0];
        $seq_id =~ s/\s//g;
        
        ## create this sequence, an analysis link, and a feature table
        $seq = $doc->createAndAddSequence($seq_id, undef, '', 'dna', 'assembly');
        $seq->addBsmlLink('analysis', '#augustus_analysis', 'input_of');
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
        
        &add_feature_group();
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
        #&add_feature('CDS', $cols[3], $cols[4], $cols[6] );
        push @cds, [ 'CDS', $cols[3], $cols[4], $cols[6] ];
    } elsif ($cols[2] eq 'start_codon') {
        &add_feature('start_codon', $cols[3], $cols[4], $cols[6] );
    } elsif ($cols[2] eq 'stop_codon') {
        &add_feature('stop_codon', $cols[3], $cols[4], $cols[6] );
    } else {
        $logger->logdie("unrecognized type: $cols[2]");
    }

}

## now handle the last set
&add_feature_group();

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
    $thing->addBsmlLink('analysis', '#augustus_analysis', 'computed_by');
    $thing->addBsmlIntervalLoc($start, $stop, $strand);

    $fg->addBsmlFeatureGroupMember( $id, $idcreator->so_used($type) );
    
    ## if type is a primary_transcript we need to add a gene too
    if ($type eq 'primary_transcript') {
        $thing = $doc->createAndAddFeature( $ft, $current_transcript_id, '', $idcreator->so_used('gene') );
        $thing->addBsmlLink('analysis', '#augustus_analysis', 'computed_by');
        $thing->addBsmlIntervalLoc($start, $stop, $strand);
        $fg->addBsmlFeatureGroupMember( $current_transcript_id, $idcreator->so_used('gene') );
    }
}

sub add_feature_group {
    if ( $extend_terminal_exon ) {
        ## add 3 bases to the last CDS, if any were found
        if (scalar @cds) {

            ## if on the reverse strand, we need to take three from column 2
            ##   if on the forward add three to column 3
            ## assumes (obviously) that all CDS in this group are on the same strand
            if ($cds[-1][3]) {
                ## here we need to sort the CDS array because the terminal one isn't 
                ##  explicitly defined and the software can write them in any order.
                ##  reverse strand, sort descending
                @cds = sort { $b->[1] <=> $a->[1] } @cds;

                $logger->debug("manually shifting 3 from reverse CDS coordinate  $cds[-1][1] on $seq_id\n") if $logger->is_debug();
                $cds[-1][1] -= 3;
            } else {
                ## here we need to sort the CDS array because the terminal one isn't 
                ##  explicitly defined and the software can write them in any order.
                ##  reverse strand, sort descending
                @cds = sort { $a->[2] <=> $b->[2] } @cds;

                $logger->debug("manually pushing 3 onto forward CDS coordinate  $cds[-1][2] on $seq_id\n") if $logger->is_debug();
                $cds[-1][2] += 3;
            }

            for my $cd ( @cds ) {
                &add_feature( @{$cd} );
            }

            undef @cds;
        }
    }

    ## pull a new gene id
    $current_transcript_id = $idcreator->new_id( db      => $options{project},
                                                 so_type => 'gene',
                                                 prefix  => $options{command_id}
                                               );
    $fg = $doc->createAndAddFeatureGroup( $seq, '', $current_transcript_id );
}

sub check_parameters {
    
    ## make sure input file exists
    if (! -e $options{'input'}) { $logger->logdie("input file $options{'input'} does not exist") }
    
    $options{'fasta_file'} = '' unless ($options{'fasta_file'});
    $options{'project'}    = 'unknown' unless ($options{'project'});
    $options{'command_id'} = '0' unless ($options{'command_id'});
    
    return 1;
}

