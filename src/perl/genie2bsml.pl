#!/usr/local/bin/perl

=head1  NAME 

genie2bsml.pl - convert genie GFF output to BSML

=head1 SYNOPSIS

USAGE: genie2bsml.pl 
        --input=/path/to/genie.output.file.gff 
        --output=/path/to/output.bsml
      [ --project=aa1 
        --fasta_file=/path/to/somefile.fsa 
      ]

=head1 OPTIONS

B<--input,-i> 
    Input file file from a genie scan.

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

This script is used to convert the output from a genie search into BSML.

=head1 INPUT

You define the input file using the --input option.  This file does not need any
special file extension.  The regular output of genie looks like this:

    aa1.assembly.19151	genie	IG	2	9465	1.602351102e+152	+	0
    aa1.assembly.19151	genie	CDS	9469	9796	464664919.7	-	0	transcript_id 1
    aa1.assembly.19151	genie	Exon	9466	9796	8730446036	-	0	transcript_id 1
    aa1.assembly.19151	genie	Stop	9469	9469	464664919.7	-	0	transcript_id 1
    aa1.assembly.19151	genie	Splice3	9797	9797	15935.71776	-	2	transcript_id 1
    aa1.assembly.19151	genie	Intron	9797	9855	1655.533863	-	2	transcript_id 1

which has the general format of:

    <seqname> <source> <feature> <start> <end> <score> <strand> <frame> [group]

This script currently assumes that the input GFF file only contains the output from
the analysis of a single FASTA sequence.  

=head1 OUTPUT


Base positions from the input file are renumbered so that positions start at zero.  Also,
the stop coordinates of each terminal CDS is extended 3bp to include the stop codon.

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
## we have to hold each cds in an array so that we can add 3 bases to the last
##  one (genie doesn't include the stop codon within the cds)
my @cds;

## go through the file
while (<$ifh>) {
    chomp;
    my @cols = split(/\t/);
    
    ## has the sequence been defined yet?  it's in the first column
    ##  this should happen on the first row only
    unless ($seq_id) {
        $seq_id = $cols[0];
        $seq_id =~ s/\s//g;
        $logger->debug("processing seq_id: $seq_id\n") if $logger->is_debug();
        
        ## create this sequence, an analysis link, and a feature table
        $seq = $doc->createAndAddSequence($seq_id);
        $seq->addBsmlLink('analysis', '#genie_analysis');
        $ft = $doc->createAndAddFeatureTable($seq);
        
        ##  also add a link to the fasta file (Seq-data-import) if requested
        if ($options{'fasta_file'}) {
            $logger->debug("adding link to fasta_file: $options{'fasta_file'}\n") if $logger->is_debug();
            $doc->createAndAddSeqDataImport( $seq, 'fasta', $options{'fasta_file'}, '', $seq_id);
        }
    }

    $current_group_name = $cols[8] || '';

    ## if column 9 (group) is defined and is different than the last one we need to
    ##  create a new feature group
    if ($current_group_name && $current_group_name ne $last_group_name) {
        
        ## remember this group name
        $last_group_name = $current_group_name;
        
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
        
        ## pull a new gene id (in genie this = primary transcript)
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
    

    ## Exon
    if ($cols[2] eq 'Exon') {
        &add_feature('exon', $cols[3], $cols[4], $cols[6], ($cols[8] || 0) );

    ## Intron
    } elsif ($cols[2] eq 'Intron') {
        &add_feature('intron', $cols[3], $cols[4], $cols[6], ($cols[8] || 0) );

    ## Splice3
    } elsif ($cols[2] eq 'Splice3') {
        &add_feature('splice_site', $cols[3], $cols[4], $cols[6], ($cols[8] || 0) );

    ## Splice5
    } elsif ($cols[2] eq 'Splice5') {
        &add_feature('splice_site', $cols[3], $cols[4], $cols[6], ($cols[8] || 0) );

    ## CDS
    } elsif ($cols[2] eq 'CDS') {
        push @cds, [ 'CDS', $cols[3], $cols[4], $cols[6], ($cols[8] || 0) ];
        #&add_feature('CDS', $cols[3], $cols[4], $cols[6], ($cols[8] || 0) );

    ## Prim_Trans
    } elsif ($cols[2] eq 'Prim_Trans') {
        &add_feature('primary_transcript', $cols[3], $cols[4], $cols[6], ($cols[8] || 0) );

    ## Stop
    } elsif ($cols[2] eq 'Stop') {
        &add_feature('transcription_end_site', $cols[3], $cols[4], $cols[6], ($cols[8] || 0) );

    ## Start
    } elsif ($cols[2] eq 'Start') {
        &add_feature('transcription_start_site', $cols[3], $cols[4], $cols[6], ($cols[8] || 0) );

    ## UTR5
    } elsif ($cols[2] eq 'UTR5') {
        &add_feature('five_prime_UTR', $cols[3], $cols[4], $cols[6], ($cols[8] || 0) );

    ## UTR3
    } elsif ($cols[2] eq 'UTR3') {
        &add_feature('three_prime_UTR', $cols[3], $cols[4], $cols[6], ($cols[8] || 0) );

    ## IG
    } elsif ($cols[2] eq 'IG') {
        &add_feature('intergenic_region', $cols[3], $cols[4], $cols[6], ($cols[8] || 0) );

    } else {
        $logger->logdie("unrecognized type: $cols[2]");
    }

}


## add the analysis element
$doc->createAndAddAnalysis(
                            id => 'genie_analysis',
                            sourcename => $options{'output'},
                          );

## now write the doc
$doc->write($options{'output'});

exit;

sub add_feature {
    my ($type, $start, $stop, $strand, $group) = @_;
    
    $logger->debug("add_feature($type, $start, $stop, $strand, $group)\n") if $logger->is_debug();
    
    $id = $idcreator->new_id( db => $options{project}, so_type => $type, prefix => $options{command_id} );
    $thing = $doc->createAndAddFeature( $ft, $id, '', $idcreator->so_used($type) );
    $thing->addBsmlLink('analysis', '#genie_analysis');
    $thing->addBsmlIntervalLoc($start, $stop, $strand);
    
    ## some features aren't added to a group
    if ($group) {
        $fg->addBsmlFeatureGroupMember( $id, $idcreator->so_used($type) );
    }
    
    ## if type is a primary_transcript we need to add a gene too
    if ($type eq 'primary_transcript') {
        $thing = $doc->createAndAddFeature( $ft, $current_transcript_id, '', $idcreator->so_used('gene') );
        $thing->addBsmlLink('analysis', '#genie_analysis');
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

