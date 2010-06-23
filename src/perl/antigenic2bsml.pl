#!/usr/bin/perl

=head1  NAME 

antigenic2bsml.pl - convert ANTIGENIC output to BSML

=head1 SYNOPSIS

USAGE: antigenic2bsml.pl --input=/path/to/antigenic_file --output=/path/to/output.bsml --project=aa1 --input=/path/to/antigenic/input.fsa --compress_bsml_output=0

=head1 OPTIONS

B<--input,-i> 
    Input file file from a antigenic run [long format only!].

B<--project>
    [OPTIONAL]
    The project (used for id making).
    Default: unknown

B<--input>
    The fasta file used as input to the antigenic run.

B<--compress_bsml_output>
    Will create gzipped output

B<--id_repository,-r>
    Required for creating feature identifiers.  Each project should have
    its own id_repository directory - use the full path to it here.  This
    is used by the IdGenerator.pm module.

B<--debug,-d> 
    Debug level.  Use a large number to turn on verbose debugging. 

B<--log,-l> 
    Log file

B<--output,-o> 
    Output BSML file (will be created, must not exist)

B<--help,-h> 
    This help message

=head1   DESCRIPTION

This script is used to convert the output from ANTIGENIC into BSML.

=head1 INPUT

ANTIGENIC can be run using multiple input sequences simultaneously, and this
script supports parsing single or multiple input result sets.
 
You define the input file using the --input option.  This file does not need any
special file extension.

If the file does not exist, the program will automatically search for a gzip'ed version 
(same filename, with the extension '.gz').  If a '.gz' version exists, will read in the 
gzip'ed file.

=head1 OUTPUT

After parsing the input file, a file specified by the --output option is created.  This script
will fail if it already exists.

=head1 CONTACT

David Riley
driley@som.igs.umaryland.edu

=cut

use strict;
use Log::Log4perl qw(get_logger);
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use BSML::BsmlBuilder;
use BSML::BsmlParserTwig;
use BSML::BsmlRepository;
use Ergatis::IdGenerator;
use Pod::Usage;
use Ergatis::Logger;

my %options = ();
my $results = GetOptions (\%options, 
              'input|i=s',
              'output|o=s',
              'debug|d=s',
              'id_repository|r=s',
              'input=s',
              'compress_bsml_output=s',
              'project|p=s',
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

my $defline;
my $fastaFile;
my $gzip = 0;

## make sure all passed options are peachy
&check_parameters(\%options);

## we want to create ids unique to this document, which will be replaced later.  they must
##  contain the prefix that will be used to look up a real id, such as ath1.gen.15
my $next_id = 1;

## we want a new doc
my $doc = new BSML::BsmlBuilder();

## we're going to generate ids
my $idgen = Ergatis::IdGenerator->new( id_repository => $options{id_repository} );

## open the input file for parsing (even if it's gziped (bug 2591))
my $mode = "<";
$mode .= ":gzip" if($options{'input'} =~ /.gz$/);    
open (my $ifh, $mode, $options{'input'}) || $logger->logdie("can't open input file for reading");

my @sequence_ids;
my %result_ref_hash;
#my $temp;
my %property_hash;
my $epitope_hash;
my %qualifier_hash;
my $skip_flag = 0;
my $in_epitope = 0;
my $start = 0;
my $end = 0;
my $strand;
my $score;
my $seq_name = '';
my $seq_name_to_feat_table;

#my $line;
while (<$ifh>) {

    if(/Name:\s(\S+)/) {
        $seq_name=$1;
    }
    elsif(/Start:\s(\S+)/) {
        $start=$1;
          }
    elsif(/End:\s(\S+)/) {
        $end = $1;
    }
    elsif(/Strand:\s(\S+)/) {
        $strand=$1;
    }
    elsif(/Score:\s(\S+)/) {
        $score=$1;
          }
    elsif(/Max_score_pos:\s(\S+)/) {

        my $seq;
        if( !( $seq = $doc->returnBsmlSequenceByIDR($seq_name)) ){
            $seq = $doc->createAndAddSequence($seq_name, $seq_name, undef, 'aa', 'polypeptide');
            $doc->createAndAddBsmlAttribute( $seq, 'defline', $defline);
#            $doc->createAndAddSeqDataImport( $seq, 'fasta', $options{input}, "sdi_" . $s . '_seq', $s );
            $seq->addBsmlLink('analysis', '#antigenic_analysis', 'input_of');
            $seq_name_to_feat_table->{$seq_name} = $doc->createAndAddFeatureTable($seq);
        }
        
        my $project = $options{'project'};
        if($seq_name =~ /^([^\.]+)\./ && ($options{'project'} eq 'unknown')) {
            $project = $1;
        }
        my $feature_table = $seq_name_to_feat_table->{$seq_name};
        my $ep = $doc->createAndAddFeature(
            $feature_table,
            $idgen->next_id( project => $project,
                             type => 'epitope'),
            'epitope');
        $ep->addBsmlLink('analysis', '#antigenic_analysis', 'computed_by');
        $ep->addBsmlIntervalLoc(
            $start,
            $end,
            $strand);
    }
}
my $analysis = $doc->createAndAddAnalysis(
    id => 'antigenic_analysis',
    sourcename => $options{'output'},
    program => 'antigenic',
    algorithm => 'antigenic',
    );
$doc->write($options{'output'},,$gzip);

exit;

sub check_parameters {

    ## required params
    my @required = qw( input output id_repository );
    for ( @required ) {
        if (! defined $options{$_}) {
            $logger->logdie( "$_ is a required option" );
        }
    }
    
    ## make sure input file exists
    if (! -e $options{'input'}) {
        unless(-e $options{'input'}.".gz") {
            $options{'input'}.=".gz";
        } else {
            pod2usage({-message => "Input file '$options{'input'}' does not exist!"});
        }
    } 
    
    ## make sure the fasta input was provided. (Bug 3781)
    if($options{'input'}) {
        $fastaFile = $options{'input'};
        open(IN, "$fastaFile") 
            or pod2usage({-message => "Could not open input $fastaFile ($!)"});
        while(<IN>) {
            chomp;
            if(/^>(.*)/) {
                $defline = $1;
                last;
            }
        }
        close(IN);
    } else {
        pod2usage({-message => "Option input is required"});
    }

    ## in case they want to compress the bsml output (bug 2591)
    if($options{'compress_bsml_output'}) {
        $gzip = 1;
    }
    
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

