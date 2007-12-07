#!/usr/bin/perl

#eval 'exec /usr/local/packages/perl-5.8.5/bin/perl  -S $0 ${1+"$@"}'
#    if 0; # not running under some shell
#BEGIN{foreach (@INC) {s/\/usr\/local\/packages/\/local\/platform/}};
use lib (@INC,$ENV{"PERL_MOD_DIR"});
#no lib "$ENV{PERL_MOD_DIR}/i686-linux";
#no lib ".";

# And remove this:
#use lib (@INC, '/usr/local/devel/ANNOTATION/ard/current/lib/site_perl/5.8.5/');

=head1  NAME 

snap2bsml.pl - convert snap ZFF output to BSML

=head1 SYNOPSIS

USAGE: snap2bsml.pl 
        --input=/path/to/snap.output.file.raw 
        --output=/path/to/output.bsml
      [ --project=aa1 
        --fasta_file=/path/to/somefile.fsa 
      ]

=head1 OPTIONS

B<--input,-i> 
    Input file from a snap search.

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

This script is used to convert the output from a snap search into BSML.

=head1 INPUT

You define the input file using the --input option.  This file does not need any
special file extension.  The regular output of snap looks like this:

        >Y73E7A.6
        Einit    201    325   +    90   0   2   1   Y73E7A.6
        Eterm   2175   2319   +   295   1   0   2   Y73E7A.6
        >Y73E7A.7
        Einit    201    462   +   263   0   1   1   Y73E7A.7
        Exon    1803   2031   +   379   2   2   0   Y73E7A.7
        Exon    2929   3031   +   236   1   0   0   Y73E7A.7
        Exon    3467   3624   +   152   0   2   0   Y73E7A.7
        Exon    4185   4406   +   225   1   2   2   Y73E7A.7
        Eterm   5103   5280   +    46   1   0   2   Y73E7A.7

which has the general format of:

><seq_name>
<Label> <Begin> <End> <strand> <score> <5'-overhang> <3'-overhang> <Frame> <Group> 

=head1 OUTPUT

Base positions from the input file are renumbered so that positions start at zero.

=head1 CONTACT

    Jason Inman
    jinman@tigr.org

=cut

use strict;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use Pod::Usage;
use Ergatis::Logger;

#use Papyrus::TempIdCreator;
use Ergatis::IdGenerator;

use BSML::GenePredictionBsml;
#use BSML::BsmlRepository;
#use BSML::BsmlBuilder;
#use BSML::BsmlParserTwig;

my %options = ();
my $results = GetOptions (\%options, 
			  'input|i=s',
                          'output|o=s',
                          'fasta_file|f=s',
                          'project|p=s',
                          'log|l=s',
                          'command_id=s',       ## passed by workflow
                          'logconf=s',          ## passed by workflow (not used)
                          'seq_class=s',
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

## go through the file
while (<$ifh>) {
    ## skip comment lines
    next if (/^\#/);
    chomp;
    
    ## has the sequence been defined yet?  it's in the first column
    ##  this should happen on the first row only
    if (/^>(.*)/) {
        $seq_id = $1;
        $seq_id =~ s/\s//g;
        
        ## create this sequence, an analysis link, and a feature table
        $logger->debug("adding seq_id $seq_id") if $logger->is_debug;
        $seq = $doc->createAndAddSequence($seq_id, undef, '', 'dna', $options{'seq_class'});
        $seq->addBsmlLink('analysis', '#snap_analysis');
        $ft = $doc->createAndAddFeatureTable($seq);
        
        ##  also add a link to the fasta file (Seq-data-import) if requested
        if ($options{'fasta_file'}) {
            $doc->createAndAddSeqDataImport( $seq, 'fasta', $options{'fasta_file'}, '', $seq_id);
        }
        next; # Get out of this loop.  We don't need to do anymore with this line.
    }

    # If we've reached this point, we're in a data line.  Let's begin parsing the columns.
    # Start with the group name.
    my @cols = split(/\t/);
    if ($cols[8] =~ /($seq_id.*)/) {
        $current_group_name = $1;
    } else {
        $logger->logdie("unrecognized format in column 9: $cols[8]");
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
        $logger->debug("adding new feature group with parent $current_transcript_id") if $logger->is_debug;
        $fg = $doc->createAndAddFeatureGroup( $seq, '', $current_transcript_id );
    }
    
    ## adjust both positions so that we are numbering from zero
    $cols[1]--;
    
    ## change the + and - symbols in strand column to 0 and 1, respectively
    if ($cols[3] eq '+') {
        $cols[3] = 0;
    } elsif ($cols[3] eq '-') {
        $cols[3] = 1;
    } else {
        $logger->logdie("unknown value ($cols[3]) in strand column.  expected + or -.");
    }

    ## handle each of the types we know about
    if ($cols[0] eq 'Einit' || $cols[0] eq 'Exon' || $cols[0] eq 'Esngl' || $cols[0] eq 'Eterm' ) {
        &add_feature('exon', $cols[1], $cols[2], $cols[3] );
        &add_feature('CDS',  $cols[1], $cols[2], $cols[3] );
    } else {
        $logger->logdie("unrecognized type: $cols[0]");
    }

}


## add the analysis element
$doc->createAndAddAnalysis(
                            id => 'snap_analysis',
                            sourcename => $options{'output'},
                          );

## now write the doc
$doc->write($options{'output'});

exit;

sub add_feature {
    my ($type, $start, $stop, $strand) = @_;
    
    $id = $idcreator->new_id( db => $options{project}, so_type => $type, prefix => $options{command_id} );
    $thing = $doc->createAndAddFeature( $ft, $id, '', $idcreator->so_used($type) );
    $thing->addBsmlLink('analysis', '#snap_analysis', 'computed_by');
    $thing->addBsmlIntervalLoc($start, $stop, $strand);

    $fg->addBsmlFeatureGroupMember( $id, $idcreator->so_used($type) );
    
    ## if type is a primary_transcript we need to add a gene too
    if ($type eq 'primary_transcript') {
        $thing = $doc->createAndAddFeature( $ft, $current_transcript_id, '', $idcreator->so_used('gene') );
        $thing->addBsmlLink('analysis', '#snap_analysis', 'computed_by');
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

