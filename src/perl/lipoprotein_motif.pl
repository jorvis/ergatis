#!/usr/bin/perl

BEGIN{foreach (@INC) {s/\/usr\/local\/packages/\/local\/platform/}};
use lib (@INC,$ENV{"PERL_MOD_DIR"});
no lib "$ENV{PERL_MOD_DIR}/i686-linux";
no lib ".";

=head1 NAME

lipoprotein_motif.pl - Scans for membrane lipoprotein lipid attachment sites on amino acid sequence.
    Uses prosite motif.

=head1 SYNOPSIS

USAGE: lipoprotein_motif.pl 
    --input=/path/to/some/input.fsa
    --output=/path/to/some/output.bsml
    --id_repository=/path/to/id_repository
    --project=aa1
  [ --gzip_output=1
    --is_mycoplasm=1
    --log=/path/to/some/log.file
    --debug=3
    --help ]

=head1 OPTIONS

B<--input,-i>
    Fasta file containing an amino acid sequence

B<--output,-o>
    Output bsml of predictions

B<--id_repository,-r>
    Valid id_repository (See Ergatis::IdGenerator for details)

B<--project,-p>
    Used for id generation, passed into IdGenerator

B<--gzip_output,-g>
    Compress bsml output

B<--is_mycoplasm,-m>
    Searches with a slightly different motif for species of the mycoplasm genus.

B<--log,-l>
    A log file that will contain information about the run.

B<--debug,-d>
    A larger number is more verbose.

B<--help,-h>
    Will print this message.

=head1  DESCRIPTION

    Using a regular expression (from $SGC_SCRIPTS/lipoprotein_update.dbi) this program will search
    amino acid sequence for the membrane lipoprotein lipid attachment site.  Only the very begining of
    the sequence is searched.

=head1  INPUT

    Input if a normal fasta file containing amino acid sequence.

=head1  OUTPUT
    
    Outputs a bsml file.
    
=head1  CONTACT

    Kevin Galens
    kgalens@tigr.org

=cut

use strict;
use warnings;
use Pod::Usage;
use Ergatis::IdGenerator;
use Ergatis::Logger;
use BSML::BsmlBuilder;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);

#####GLOBALS#########
use constant FEAT_TITLE => "PDOC00013 :: Prokaryotic membrane lipoprotein lipid attachment site";
my $inFileName;
my $outFileName;
my $regex;
my $mycoplasm = 0;
my $doc = new BSML::BsmlBuilder();
my %seqsAdded;
my $idMaker;
my $gzip;
my $project;
#####################

my %options = ();
my $results = GetOptions (\%options, 
                          'input|i=s',
                          'output|o=s',
                          'is_mycoplasm|m=s',
                          'id_repository|r=s',
                          'project|p=s',
                          'gzip_output|g=s',
                          'debug|d=s',
                          'log|l=s',
                          'help|h') || &_pod;

#Setup the logger
my $logfile = $options{'log'} || Ergatis::Logger::get_default_logfilename();
my $logger = new Ergatis::Logger('LOG_FILE'=>$logfile,
				  'LOG_LEVEL'=>$options{'debug'});
$logger = $logger->get_logger();

#Make sure the input options were okay
&checkOptions(\%options);

#Set the actual motif (regex).
$regex = &setMotif($mycoplasm);

#Get the sequence from the fasta file.
my $sequence = &getSequence($inFileName);

#Check for matches
foreach my $defline (keys %{$sequence}) {

    #Each protein sequence is stored in the sequence hash with the defline of the fasta entry as the
    #key.  This is incase there are multiple sequences in the fasta file. 
    my $protein = $sequence->{$defline};

    #Actually use the regex and make sure that the motif is long enough, but didn't occur to
    #far in the sequence.
    if($protein =~ /($regex)/g && (pos $protein > 14 && pos $protein < 36) ) {

        #Store the hit in bsml.
        &bsml_storeHit($defline, $protein, $1);
    }
}

#Add the analysis
$doc->createAndAddAnalysis( 'id' => 'lipoprotein_motif_analysis',
                            'sourcename' => $inFileName );

#Write the bsml file.
$doc->write($outFileName, '', $gzip);

#####################SUB ROUTINES#################################################

#Description:  Will take in a definition line, protein sequence, and the signal
#              that was found to match.  It will then find the coordinates of that
#              match and store it as a feature (of type signal_peptide).  
#Parameters:   $defline - definition line for input fasta sequence
#              $protein - string contiaining protein sequence
#              $match - the signal that was found.
#Returns:      Nothing
sub bsml_storeHit {
    my ($defline, $protein, $match) = @_;

    my $id;
    $id = $1 if($defline =~ /^([^\s]+)/);
    &_die("Couldn't parse out the id out of the header (>$defline)") unless($id);

    $project = "";
    $project = $1 if($id =~ /^([^\.]+)\./);
    &_die("Could not parse project out of id '$id'") unless($project);
    
    my $seq = $doc->returnBsmlSequenceByIDR( $id );

    unless($seq) {
        $seq = $doc->createAndAddSequence( $id, $id, length($protein), 'aa', 'polypeptide' );
        $doc->createAndAddBsmlAttribute( $seq, 'defline', $defline);
        $doc->createAndAddSeqDataImport( $seq, 'fasta', $inFileName, '', $id ); 
        $doc->createAndAddLink( $seq, 'analysis', '#lipoprotein_motif_analysis', '#input_of' );
    }

    my $fTable = $doc->createAndAddFeatureTable( $seq );
    my $featId = $idMaker->next_id('type'    => 'signal_peptide',
                                  'project' => $project );

    my $feat = $doc->createAndAddFeature( $fTable, $featId, FEAT_TITLE, 'signal_peptide');

    my $start = index($protein, $match);
    my $end = $start + length($match);
    my $comp;
    ($start, $end, $comp) = ($start > $end) ? ($end, $start, 1) : ($start, $end, 0);

    $doc->createAndAddIntervalLoc( $feat, $start, $end, $comp);
    $logger->logdie("Feat was not set") unless($feat);
    $doc->createAndAddLink( $feat, 'analysis', '#lipoprotein_motif_anlaysis', 'computed_by');

}

#Description:  Will retrieve sequences from a fasta file.  
#Parameters:   Input fasta file name.
#Returns:      Hash reference - containing the fasta definition lines as the key
#                     and protein sequence as the value.
sub getSequence {
    my $input = shift;
    my $tmpSeq = "";
    my $retSeq = {};
    my $id;
    open(IN, "< $input") or
        &_die("Unable to open input file $input ($!)");
    
    while(<IN>) {
        chomp;
        if(/^>(.*)/) {
            $retSeq->{$id} =  $tmpSeq if($tmpSeq);
            $id = $1;
            $tmpSeq = "";
        } elsif($_ !~ /^\s+$/) {
            $tmpSeq.=$_;
        }
    }
    $retSeq->{$id} = $tmpSeq if($tmpSeq);

    return $retSeq;
    
}

#Description:  Chooses which regular expression to use for the motif.
#Parameters:   Boolean - A non zero value indicates a species of the Mycoplasm genus.
#Returns:      The regular expression to be used to find the lipoprotein signal.
sub setMotif {
    my $myco = shift;
    my $motif = ($myco) ? 
        "((^.{0,6}[KR]).{0,18}[^DERK][^DERK][^DERK][^DERK][^DERK][^DERK][LIVMFWSTAG][LIVMFWSTAG]".
        "[LIVMFWSTAGCQY][AGSTKRQ]C)" :
        "((^.{0,6}[KR]).{0,18}[^DERK][^DERK][^DERK][^DERK][^DERK][^DERK][LIVMFWSTAG][LIVMFWSTAG]".
        "[LIVMFWSTAGCQY][AGS]C)";
    return $motif;
}

#This checks the options.
sub checkOptions {
    my $opt = shift;
    my $err = "";

    if(!$opt->{'input'}) {
        $err.= "Input option is required\n";
    } else {
        $inFileName = $opt->{'input'};
        $err.="Input file $inFileName does not exist\n" unless(-e $inFileName);
    }

    if(!$opt->{'output'}) {
        $err.="Output option is required\n";
    } else {
        $outFileName = $opt->{'output'};
    }

    if($opt->{'is_mycoplasm'}) {
        $mycoplasm = 1;
    }

    if($opt->{'id_repository'}) {
        $idMaker = new Ergatis::IdGenerator( 'id_repository' => $opt->{'id_repository'} );
        $idMaker->set_pool_size( 'signal_peptide' => 25 );
        $project = $opt->{'project'};
    } else {
        $err.= "id_repository option is required\n" unless($opt->{'id_repository'});
    }

    if($opt->{'gzip_output'}) {
        $gzip = 1;
    } else {
        $gzip = 0;
    }


    &_die($err) if($err);
}

# Pod.  Saves some typing.
sub _pod {
    pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDERR} );
}

#Die.  Make sure I use logger everytime.
sub _die {
    my $msg = shift;
    $logger->logdie($msg);
}
##EOF
