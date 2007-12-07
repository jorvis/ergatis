#!/usr/bin/perl
BEGIN{foreach (@INC) {s/\/usr\/local\/packages/\/local\/platform/}};
use lib (@INC,$ENV{"PERL_MOD_DIR"});
no lib "$ENV{PERL_MOD_DIR}/i686-linux";
no lib ".";

=head1  NAME 

phat2bsml.pl - convert phat GFF output to BSML

=head1 SYNOPSIS

USAGE: phat2bsml.pl 
        --input=/path/to/phat.output.file.gff 
        --output=/path/to/output.bsml
      [ --project=aa1 ]

=head1 OPTIONS

B<--input,-i> 
    Input file from a phat scan.

B<--debug,-d> 
    Debug level.  Use a large number to turn on verbose debugging. 

B<--project,-p> 
    Project ID.  Used in creating feature ids.  Defaults to 'unknown' if
    not passed.

B<--command_id> 
    This is passed automatically by workflow and is used in creating unique
    feature IDs.  It represents the command ID of the phat2bsml step withing the
    phat workflow component.

B<--log,-l> 
    Log file

B<--output,-o> 
    Output BSML file (will be created, must not exist)

B<--help,-h> 
    This help message

=head1   DESCRIPTION

This script is used to convert the output from a phat search into BSML.

=head1 INPUT

You define the input file using the --input option.  This file does not need any
special file extension.  The regular output of phat looks like this:

    sma1.assembly.29055     FullPhat        CDS     8194    8426    1.00    -       1       gene "Phat1" frame"1"
    sma1.assembly.29055     FullPhat        CDS     10940   11067   0.35    -       0       gene "Phat1" frame"0"
    sma1.assembly.29055     FullPhat        CDS     26248   26370   0.52    +       0       gene "Phat2" frame"0"
    sma1.assembly.29055     FullPhat        CDS     27212   27390   1.00    +       0       gene "Phat2" frame"1"
    sma1.assembly.29055     FullPhat        CDS     29317   29539   0.99    +       1       gene "Phat2" frame"1"
    sma1.assembly.29055     FullPhat        CDS     30204   30371   0.79    +       0       gene "Phat2" frame"2"
    sma1.assembly.29055     FullPhat        CDS     33893   34008   0.60    +       0       gene "Phat3" frame"1"
    sma1.assembly.29055     FullPhat        CDS     35857   36031   0.97    +       1       gene "Phat3" frame"1"
    sma1.assembly.29055     FullPhat        CDS     47539   47666   0.75    +       0       gene "Phat3" frame"0"
    sma1.assembly.29055     FullPhat        CDS     51132   51354   1.00    +       1       gene "Phat3" frame"0"
    sma1.assembly.29055     FullPhat        CDS     54721   54870   0.97    +       0       gene "Phat3" frame"0"


which has the general format of:

    <seqname> <source> <feature> <start> <end> <score> <strand> <frame> <group>

This script currently assumes that the input GFF file only contains the output from
the analysis of a single FASTA sequence.   The raw output files from phat are not currently
used.

=head1 OUTPUT

Base positions from the input file are renumbered so that positions start at zero.  Terminal
CDS predictions from phat contain the stop codon, so no adjustment to include it within the
BSML output is made.

=head1 CONTACT

    Joshua Orvis
    jorvis@tigr.org

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
              'project|p=s',
              'log|l=s',
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
my $idcreator = new Papyrus::TempIdCreator();

## open the input file for parsing
open (my $ifh, $options{'input'}) || $logger->logdie("can't open input file for reading");

my $seq_id;
my ($seq, $ft, $fg);
my ($last_group_name, $current_group_name, $current_transcript_id);
my ($thing, $id);
## we have to hold each cds in an array so that we can add 3 bases to the last
##  one (phat doesn't include the stop codon within the cds)
my @cds;

## go through the file
while (<$ifh>) {
    chomp;
    ## skip whitespace
    next if (/^\s*$/);
    my @cols = split(/\t/);
    
    ## die if we didn't get enough columns
    unless ( scalar @cols == 9 ) {
        $logger->logdie("failed to split line $. into 9 fields.  quitting.");
    }
    
    ## has the sequence been defined yet?  it's in the first column
    ##  this should happen on the first row only
    unless ($seq_id) {
        $seq_id = $cols[0];
        $seq_id =~ s/\s//g;
        $logger->debug("processing seq_id: $seq_id\n") if $logger->is_debug();
        
        ## create this sequence, an analysis link, and a feature table
        $seq = $doc->createAndAddSequence($seq_id, undef, '', 'dna', 'assembly');
        $seq->addBsmlLink('analysis', '#phat_analysis', 'input_of');
        $ft = $doc->createAndAddFeatureTable($seq);
    }

    if ( $cols[8] =~ /gene \"(.+?)\" frame/ ) {
        $current_group_name = $1;
    } else {
        $logger->logdie("unrecognized column 9 ($cols[8]), line $. .  expected something like: gene \"Phat2\" frame\"2\"\n");
    }

    ## if a group is defined and is different than the last one we need to
    ##  create a new feature group
    if ($current_group_name && $current_group_name ne $last_group_name) {
        
        ## remember this group name
        $last_group_name = $current_group_name;
        
#        ## add 3 bases to the last CDS, if any were found
#        if (scalar @cds) {
#            
#            ## if on the reverse strand, we need to take three from column 2
#            ##   if on the forward add three to column 3
#            ## assumes (obviously) that all CDS in this group are on the same strand
#            if ($cds[-1][3]) {
#                ## here we need to sort the CDS array because the terminal one isn't 
#                ##  explicitly defined and the software can write them in any order.
#                ##  reverse strand, sort descending
#                @cds = sort { $b->[1] <=> $a->[1] } @cds;
#            
#                $logger->debug("manually shifting 3 from reverse CDS coordinate  $cds[-1][1] on $seq_id\n") if $logger->is_debug();
#                $cds[-1][1] -= 3;
#            } else {
#                ## here we need to sort the CDS array because the terminal one isn't 
#                ##  explicitly defined and the software can write them in any order.
#                ##  reverse strand, sort descending
#                @cds = sort { $a->[2] <=> $b->[2] } @cds;
#
#                $logger->debug("manually pushing 3 onto forward CDS coordinate  $cds[-1][2] on $seq_id\n") if $logger->is_debug();
#                $cds[-1][2] += 3;
#            }
            
            for my $cd ( @cds ) {
                &add_feature( @{$cd} );
            }
            
            undef @cds;
#        }
        
        ## pull a new gene id (in phat this = primary transcript)
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

    ## phat only predicts CDS features
    if ($cols[2] eq 'CDS') {
        push @cds, [ 'CDS', $cols[3], $cols[4], $cols[6], ($cols[8] || 0) ];
    } else {
        $logger->logdie("unrecognized feature type in raw phat output. expected only CDS: $cols[2]");
    }

}


## add the analysis element
$doc->createAndAddAnalysis(
                            id => 'phat_analysis',
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
    $thing->addBsmlLink('analysis', '#phat_analysis', 'computed_by');
    $thing->addBsmlIntervalLoc($start, $stop, $strand);
    
    ## some features aren't added to a group
    if ($group) {
        $fg->addBsmlFeatureGroupMember( $id, $idcreator->so_used($type) );
    }
}

sub check_parameters {
    
    ## make sure input file exists
    if (! -e $options{'input'}) { $logger->logdie("input file $options{'input'} does not exist") }
    
    $options{'project'}    = 'unknown' unless ($options{'project'});
    $options{'command_id'} = '0' unless ($options{'command_id'});
    
    return 1;
}

