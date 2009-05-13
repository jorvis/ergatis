#!/usr/bin/perl

BEGIN{foreach (@INC) {s/\/usr\/local\/packages/\/local\/platform/}};
use lib (@INC,$ENV{"PERL_MOD_DIR"});
no lib "$ENV{PERL_MOD_DIR}/i686-linux";
no lib ".";

=head1  NAME 

signalp2bsml.pl - convert SignalP output to BSML

=head1 SYNOPSIS

USAGE: signalp2bsml.pl --input=/path/to/signalp_file --output=/path/to/output.bsml

=head1 OPTIONS

B<--input,-i> 
    Input file file from a signalp run.

B<--debug,-d> 
    Debug level.  Use a large number to turn on verbose debugging. 

B<--log,-l> 
    Log file

B<--output,-o> 
    Output BSML file (will be created, must not exist)

B<--help,-h> 
    This help message

=head1   DESCRIPTION

This script is used to convert the output from SignalP into BSML.

=head1 INPUT

SignalP can be run using multiple input sequences simultaneously, and this
script supports parsing single or multiple input result sets.
Output file types from 'nn', 'hmm', or 'nn+hmm' prediction methods are all supported 
and recognized automatically. 
 
You define the input file using the --input option.  This file does not need any
special file extension.

=head1 OUTPUT

After parsing the input file, a file specified by the --output option is created.  This script
will fail if it already exists.

=head1 CONTACT

Brett Whitty
bwhitty@tigr.org

=cut

use warnings;
use strict;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use Pod::Usage;
use Ergatis::Logger;
use BSML::BsmlRepository;
use Ergatis::IdGenerator;
use BSML::BsmlBuilder;
use BSML::BsmlParserTwig;
use Data::Dumper;

my %options = ();
my $results = GetOptions (\%options, 
              'input|i=s',
              'query_file_path|q=s',
              'output|o=s',
              'debug|d=s',
              'command_id=s',       ## passed by workflow
              'logconf=s',          ## passed by workflow (not used)
              'project|p=s',
              'id_repository=s',
              'log|l=s',
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

## we want to create ids unique to this document, which will be replaced later.  they must
##  contain the prefix that will be used to look up a real id, such as ath1.gen.15
my $next_id = 1;

## we want a new doc
my $doc = new BSML::BsmlBuilder();

## we're going to generate ids
my $idcreator = new Ergatis::IdGenerator( id_repository => $options{'id_repository'});
$idcreator->set_pool_size( 'signal_peptide' => 25 );

## open the input file for parsing
open (my $ifh, $options{'input'}) || $logger->logdie("can't open input file for reading");

my @sequence_ids;
my %deflines;
my %result_ref_hash;
my $temp;
my $next_seq_flag = 0;
while (<$ifh>) {
    chomp;

    ##recognize start of results for an individual sequence
    if ($next_seq_flag || /^-{70}/) {
    $logger->debug("start of sequence results entry found") if($logger->is_debug());  
    my %nn_results;
    my %hmm_results;
    my $defline;
    if ($next_seq_flag) {
        $next_seq_flag = 0;
        $defline = $_;
    } else {
        $defline = <$ifh>;
    }
    chomp $defline;
    $defline =~ s/^>//;
    $defline =~ /^(\S+)/;
    my $sequence_id = $1;
    $logger->debug("sequence id is '$sequence_id'") if($logger->is_debug());  
    push(@sequence_ids, $sequence_id);
    $deflines{$sequence_id} = $defline;    
    
    while (my $result_line = <$ifh>) {
        chomp $result_line;
        if ($result_line =~ /^SignalP-NN result:/) {
            $logger->debug("parsing signalp nn results") if($logger->is_debug());  
            $temp = <$ifh>;
            chomp $temp;
            $temp =~ /length = (\d+)/ || $logger->error("Couldn't parse sequence length in NN result.") if ($logger->is_error);
            $temp = <$ifh>; # we will discard the header line
            $temp = <$ifh>;
            chomp $temp;
            $temp =~ /^  max\. C\s+([^\s]+)\s+([^\s]+)\s+([^\s]+)\s+([^\s]+)/ || $logger->logdie("Couldn't parse max. C line in NN result");
            $nn_results{'maxc_pos'} = $1;
            $nn_results{'maxc_value'} = $2;
            $nn_results{'maxc_cutoff'} = $3;
            $nn_results{'maxc_signal_peptide'} = $4;
            
            $temp = <$ifh>;
            chomp $temp;
            $temp =~ /^  max\. Y\s+([^\s]+)\s+([^\s]+)\s+([^\s]+)\s+([^\s]+)/ || $logger->logdie("Couldn't parse max. Y line in NN result");
            $nn_results{'maxy_pos'} = $1;
            $nn_results{'maxy_value'} = $2;
            $nn_results{'maxy_cutoff'} = $3;
            $nn_results{'maxy_signal_peptide'} = $4;
            
            $temp = <$ifh>;
            chomp $temp;
            $temp =~ /^  max\. S\s+([^\s]+)\s+([^\s]+)\s+([^\s]+)\s+([^\s]+)/ || $logger->logdie("Couldn't parse max. S line in NN result");
            $nn_results{'maxs_pos'} = $1;
            $nn_results{'maxs_value'} = $2;
            $nn_results{'maxs_cutoff'} = $3;
            $nn_results{'maxs_signal_peptide'} = $4;
            
            $temp = <$ifh>;
            chomp $temp;
            $temp =~ /^  mean S\s+([^\s]+)\s+([^\s]+)\s+([^\s]+)\s+([^\s]+)/ || $logger->logdie("Couldn't parse mean S line in NN result");
            $nn_results{'means_pos'} = $1;
            $nn_results{'means_value'} = $2;
            $nn_results{'means_cutoff'} = $3;
            $nn_results{'means_signal_peptide'} = $4;
            
            $temp = <$ifh>;
            chomp $temp;
            $temp =~ /^       D\s+([^\s]+)\s+([^\s]+)\s+([^\s]+)\s+([^\s]+)/ || $logger->logdie("Couldn't parse D line in NN result");
            $nn_results{'d_pos'} = $1;
            $nn_results{'d_value'} = $2;
            $nn_results{'d_cutoff'} = $3;
            $nn_results{'d_signal_peptide'} = $4;
            
            $temp = <$ifh>;
            chomp $temp;
            if ($temp =~ /^# Most likely cleavage site between pos\. (\d+) and (\d+): ([A-Z\-]+)$/) {
                $nn_results{'cleavage_site_position'} = $1."-".$2;
            }
            $result_ref_hash{$sequence_id, 'nn'} = \%nn_results;

        } elsif ($result_line =~ /^SignalP-HMM result:/) {
            $logger->debug("parsing signalp hmm results") if($logger->is_debug());  
            $temp = <$ifh>; # we will discard the header line
            $temp = <$ifh>;
            chomp $temp;
            ## If HMM predicts no signal peptide, prediction will be 'Non-secretory protein'
            $temp =~ /^Prediction: (.*)$/ || $logger->logdie("Couldn't parse prediction line in HMM result");
            $hmm_results{'prediction'} = $1;
            $temp = <$ifh>;
            chomp $temp;
            $temp =~ /^Signal peptide probability: (.*)$/ || $logger->logdie("Couldn't parse signal peptide probability line in HMM result");
            $hmm_results{'signal_peptide_probability'} = $1;
            $temp = <$ifh>;
            chomp $temp;
            if ($temp =~ /^Signal anchor probability: (.*)$/) {
                $hmm_results{'signal_anchor_probability'} = $1;
                $temp = <$ifh>;
                chomp $temp;            
            }
            $temp =~ /^Max cleavage site probability: ([0-9]\.[0-9]{3}) between pos\. ([0-9\-]+)\s+and\s+([0-9\-]+)/ || $logger->logdie("Couldn't parse max cleavage site probability line in HMM result");
            $hmm_results{'max_cleavage_site_probability'} = $1;
            $hmm_results{'cleavage_site_position'} = $2."-".$3;

            $result_ref_hash{$sequence_id, 'hmm'} = \%hmm_results;
            last;
        } elsif ($result_line =~ /^-{70}/) { ## we'll hit this if result file contains only nn results
            $logger->debug("found beginning of new results block after nn results\nshould only occur if run method was 'nn' only") if($logger->is_debug());  
            $next_seq_flag = 1; ## so flag that we've hit the start of next sequence's results
            last;
        }
    }
    }
}

## map these hash keys to ontology terms
my %nn_onto = (
                'maxs_pos'              =>  's-position',
                'maxs_value'            =>  's-score',
                'maxs_cutoff'           =>  's-cutoff',
                'maxs_signal_peptide'   =>  's-prediction',
                'means_pos'             =>  's-mean-position',
                'means_value'           =>  's-mean',
                'means_cutoff'          =>  's-mean-cutoff',
                'means_signal_peptide'  =>  's-mean-prediction',
                'd_pos'                 =>  'd-position',
                'd_value'               =>  'd-score',
                'd_cutoff'              =>  'd-cutoff',
                'd_signal_peptide'      =>  'd-prediction',
                'maxc_pos'              =>  'c-position',
                'maxc_value'            =>  'c-score',
                'maxc_cutoff'           =>  'c-cutoff',
                'maxc_signal_peptide'   =>  'c-prediction',
                'maxy_pos'              =>  'y-position',
                'maxy_value'            =>  'y-score',
                'maxy_cutoff'           =>  'y-cutoff',
                'maxy_signal_peptide'   =>  'y-prediction',
               ); 

foreach my $s(@sequence_ids) {
    my $seq = $doc->createAndAddSequence(
                                        $s,             # id
                                        $s,             # title
                                        '',             # length
                                        'aa',           # molecule
                                        'polypeptide'   # class
                                        );

    $doc->createAndAddSeqDataImport(
                                        $seq, 
                                        'fasta', 
                                        $options{'query_file_path'}, 
                                        '', 
                                        $s
                                   );

    $doc->createAndAddBsmlAttribute( 
                                        $seq, 
                                        'defline', 
                                        $deflines{$s},
                                   );

    foreach my $a(('nn','hmm')) {
        $seq->addBsmlLink('analysis', '#signalp_'.$a.'_analysis', 'input_of');
    }
    
    my $ft;

    ## create BSML for nn predictions
    if ($result_ref_hash{$s,'nn'}->{'cleavage_site_position'}) {
        if (!$ft){
            $ft  = $doc->createAndAddFeatureTable($seq);
        }
        ## parse out cleavage site coords
        my @coords = parse_position($result_ref_hash{$s,'nn'}->{'cleavage_site_position'});
        my $site_pos = shift(@coords);
        
        ## create signal peptide feature
        my $signalp_id = $idcreator->next_id( project   => $options{project},
                                              type      => 'signal_peptide',
                                            );
        my $signalp = $doc->createAndAddFeature($ft, $signalp_id, '', 'signal_peptide');
        $signalp->addBsmlIntervalLoc('0', $site_pos);
        $signalp->addBsmlLink('analysis', '#signalp_nn_analysis', 'computed_by');
        foreach my $att((
                            'maxs_pos',
                            'maxs_value',
                            'maxs_cutoff',
                            'maxs_signal_peptide',
                            'means_pos',
                            'means_value',
                            'means_cutoff',
                            'means_signal_peptide',
                            'd_pos',
                            'd_value',
                            'd_cutoff',
                            'd_signal_peptide',
                        )) {
            $doc->createAndAddBsmlAttribute(
                                        $signalp,
                                        $nn_onto{$att}, 
                                        $result_ref_hash{$s,'nn'}->{$att}
                                       );
        }
        
        ## create cleavage site feature
        my $csite_id = $idcreator->next_id( project => $options{project},
                                         type       => 'cleavage_site',
                                       );

        my $cleavage_site = $doc->createAndAddFeature($ft, $csite_id, '', 'cleavage_site');
        $cleavage_site->addBsmlLink('analysis', '#signalp_nn_analysis', 'computed_by');
        $cleavage_site->addBsmlSiteLoc($site_pos);
        $doc->createAndAddBsmlAttribute(
                                        $cleavage_site, 
                                        'max_cleavage_site_probability', 
                                        $result_ref_hash{$s,'nn'}->{'max_cleavage_site_probability'}
                                       );
        foreach my $att((
                            'maxc_pos',
                            'maxc_value',
                            'maxc_cutoff',
                            'maxc_signal_peptide',
                            'maxy_pos',
                            'maxy_value',
                            'maxy_cutoff',
                            'maxy_signal_peptide',
                        )) {
            $doc->createAndAddBsmlAttribute(
                                        $cleavage_site,
                                        $nn_onto{$att}, 
                                        $result_ref_hash{$s,'nn'}->{$att}
                                       );
        }
    }
    
    ## create BSML for hmm predictions
    if ($result_ref_hash{$s,'hmm'}->{'prediction'} ne 'Non-secretory protein') {
        if (!$ft){
            $ft  = $doc->createAndAddFeatureTable($seq);
        }
        
        
        ## parse out cleavage site coords
        my @coords = parse_position($result_ref_hash{$s,'hmm'}->{'cleavage_site_position'});
        my $site_pos = shift(@coords);
            
        ## create signal peptide feature
        my $signalp_id = $idcreator->next_id(project => $options{project}, type => 'signal_peptide');
        my $signalp = $doc->createAndAddFeature(
                                                $ft, 
                                                $signalp_id, 
                                                '', 
                                                'signal_peptide',
                                               );
        $signalp->addBsmlIntervalLoc('0', $site_pos);
        $signalp->addBsmlLink('analysis', '#signalp_hmm_analysis', 'computed_by');
        $doc->createAndAddBsmlAttribute(
                                        $signalp, 
                                        'signal_anchor', 
                                        $result_ref_hash{$s,'hmm'}->{'signal_anchor_probability'},
                                       );
        $doc->createAndAddBsmlAttribute(
                                        $signalp, 
                                        'signal_probability', 
                                        $result_ref_hash{$s,'hmm'}->{'signal_peptide_probability'}
                                       );
        $doc->createAndAddBsmlAttribute(
                                        $signalp, 
                                        'prediction', 
                                        $result_ref_hash{$s,'hmm'}->{'prediction'}
                                       );

        ## prediction can be 'Signal anchor' or 'Signal peptide', signal anchor is uncleaved signal peptide
        if ($result_ref_hash{$s,'hmm'}->{'prediction'} eq 'Signal peptide') {
            ## create cleavage site feature
            my $csite_id = $idcreator->next_id( project => $options{project},
                                             type       => 'cleavage_site',
                                           );

            my $cleavage_site = $doc->createAndAddFeature($ft, $csite_id, '', 'cleavage_site');
            $cleavage_site->addBsmlLink('analysis', '#signalp_hmm_analysis', 'computed_by');
            $cleavage_site->addBsmlSiteLoc($site_pos);
            $doc->createAndAddBsmlAttribute(
                                $cleavage_site, 
                                'max_cleavage_site_probability', 
                                $result_ref_hash{$s,'hmm'}->{'max_cleavage_site_probability'}
                                           );
        }
    }
}
    
## add the analysis element
foreach my $a(('nn','hmm')) {
    my $analysis = $doc->createAndAddAnalysis(
                            id => 'signalp_'.$a.'_analysis',
                            sourcename => $options{'output'},
                            program => "signalp_$a",
                            algorithm => "signalp_$a",
                          );
}

## now write the doc
$doc->write($options{'output'});

exit();

sub check_parameters {
    
    ## make sure input file exists
    if (! -e $options{'input'}) { $logger->logdie("input file $options{'input'} does not exist") }
   
    ##
    if (! $options{'query_file_path'}) { $logger->logdie("must provide path to query fasta file with --query_file_path") }
    
    $options{'project'}    = 'unknown' unless ($options{'project'});
    $options{'command_id'} = '0' unless ($options{'command_id'});
    
    return 1;
}

## parse_position parses the position field to give an array containing
## start [and stop] sites
sub parse_position {
    my ($position) = @_;
    if ($position =~ /^\s*([^\s]+)\s*$/) {
        return split("-", $1);
    } else {
        return ();
    }
}

