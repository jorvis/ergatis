#!/usr/bin/perl
BEGIN{foreach (@INC) {s/\/usr\/local\/packages/\/local\/platform/}};
use lib (@INC,$ENV{"PERL_MOD_DIR"});
no lib "$ENV{PERL_MOD_DIR}/i686-linux";
no lib ".";

=head1  NAME 

ngap.pl - find N-gap positions in nucleotide sequences

=head1 SYNOPSIS

USAGE: ngap.pl 
        --input=/path/to/fasta_file.fsa
        --output=/path/to/output.bsml
        --id_repository=/path/to/some/repository
      [ --project=aa1 ]

=head1 OPTIONS

B<--id_repository,-r> 
    Path to the project's ID repository, for use by the IdGenerator module. 

B<--input,-i> 
    Input FASTA or multi-FASTA format file.

B<--debug,-d> 
    Debug level.  Use a large number to turn on verbose debugging. 

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

This script is used to identify N-gap regions in FASTA format nucleotide sequences.

=head1 INPUT

The input should be an individual or multiple nucleotide sequences in one FASTA format file.
Sequence identifiers are recognized as the first string of non-whitespace characters occurring 
after the '>' character.

=head1 OUTPUT

The output of this script is a BSML format file. Each input sequence is represented as a sequence
stub with associated features of type 'gap' with interval-loc elements specifying the position
of the N-gaps. The filename used will be that specified by the --output option.  This script
will fail if the specified output file already exists.  The file is created, and temporary IDs 
are assigned for each result element.  They are only unique to the document, and will need to be 
replaced before any database insertion.

Gap positions are specified using interbase numbering.

=head1 CONTACT

Brett Whitty
bwhitty@tigr.org

=cut

use strict;
use warnings;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use Pod::Usage;
use Ergatis::Logger;
use Ergatis::IdGenerator;
use BSML::BsmlRepository;
use BSML::BsmlBuilder;
use BSML::BsmlParserTwig;

my %options = ();
my $results = GetOptions (\%options, 
              'input|i=s',
              'output|o=s',
              'id_repository|r=s',
              'project|p=s',
              'log|l=s',
              'debug=s',
              'help|h') || pod2usage();

if (scalar keys(%options) < 1) {
    pod2usage( {-exitval=>0, -verbose => 2, -output => \*STDOUT} );
}

## assign a value to project if it wasn't provided so that tempidcreator won't complain
if (!defined($options{'project'})) {
    $options{'project'} = 'null';   
}
## assign a value to command_id if it wasn't provided so that tempidcreator won't complain
if (!defined($options{'command_id'})) {
    $options{'command_id'} = '0';
}

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

## we want to creating ids unique to this document, which will be replaced later.  they must
##  contain the prefix that will be used to look up a real id, such as ath1.gen.15
my $next_id = 1;

## we want a new doc
my $doc = new BSML::BsmlBuilder();

## we're going to generate ids
my $idcreator = new Ergatis::IdGenerator( id_repository => $options{id_repository} );

## open the input file for parsing
open (IN, $options{'input'}) || $logger->logdie("can't open input file for reading");

my $seq = '';
my $seq_id = '';
my %ngaps = ();

## Read in FASTA sequences and find Ngaps
while (<IN>) {
    chomp;
    if (/^>([^\s]+)/) {
        if ($seq ne '') {
            $ngaps{$seq_id} = &find_ngaps(\$seq);
        }
        $seq = '';
        $seq_id = $1;
    } else {
        $seq .= $_;
    }
}
if ($seq ne '') {
    $ngaps{$seq_id} = &find_ngaps(\$seq);
}
close IN;

## Write BSML output file

foreach $seq_id(keys(%ngaps)) {
    my $seq_stub = $doc->createAndAddSequence(
                                                $seq_id, 
                                                undef, 
                                                undef, 
                                                'na', 
                                                'assembly'
                                             );
    $seq_stub->addBsmlLink(
                            'analysis', 
                            '#ngap_analysis', 
                            'input_of'
                          );

    my $feature_table;
    
    ## for more efficient ID generation
    if ( scalar @{$ngaps{$seq_id}} ) {
        $idcreator->set_pool_size( gap => scalar @{$ngaps{$seq_id}} );
    }
    
    foreach my $ngap_ref(@{$ngaps{$seq_id}}) {
        unless ($feature_table) {
            $feature_table = $doc->createAndAddFeatureTable($seq_stub);
        }
        my $ngap = $doc->createAndAddFeature(
                                                $feature_table, 
                                                $idcreator->next_id(
                                                                    project => $options{'project'},
                                                                    type => 'gap'
                                                                  ),
                                                '',
                                                'gap'
                                            );
                                          
        $ngap->addBsmlLink('analysis', '#ngap_analysis', 'computed_by');
        $ngap->addBsmlIntervalLoc(
                                    $ngap_ref->[0],
                                    $ngap_ref->[1],
                                    0
                                 );
    }
}

## add analysis element
my $analysis = $doc->createAndAddAnalysis(
                                            id => 'ngap_analysis',
                                            sourcename => $options{'output'},
                                         );

## write the doc
$doc->write($options{'output'});

exit(0);

sub find_ngaps {
    my ($seq_ref) = shift @_;
    my @ngaps = ();
    
    while (${$seq_ref} =~ /[N]{1,}/ig) {
        push (@ngaps, [$-[0],$+[0]]);
    }
    return \@ngaps;
}

sub check_parameters {
    ## check required options
    my @required = qw( input output id_repository );
    for ( @required ) {
        unless ( defined $options{$_} ) {
                $logger->logdie("--$_ is a required parameter");
        }
    }

    ## make sure input file exists
    if (! -e $options{'input'}) { $logger->logdie("input file $options{'input'} does not exist") }
    
    return 1;
}
