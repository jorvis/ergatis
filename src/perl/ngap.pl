#!/usr/bin/perl
use lib (@INC,$ENV{"PERL_MOD_DIR"});
no lib "$ENV{PERL_MOD_DIR}/i686-linux";
no lib ".";

=head1  NAME 

ngap.pl - find N-gap positions in nucleotide sequences

=head1 SYNOPSIS

USAGE: ngap.pl 
        --input=/path/to/fasta_file.fsa
        --output=/path/to/output.bsml
      [ --project=aa1 ]

=head1 OPTIONS

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
BEGIN {
    use Workflow::Logger;
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

## we want to creating ids unique to this document, which will be replaced later.  they must
##  contain the prefix that will be used to look up a real id, such as ath1.gen.15
my $next_id = 1;

## we want a new doc
my $doc = new BSML::BsmlBuilder();

## we're going to generate ids
my $idcreator = new Papyrus::TempIdCreator();

## open the input file for parsing
open (IN, $options{'input'}) || $logger->logdie("can't open input file for reading");

my $seq = '';
my $seq_id = '';
my %ngaps = ();

## Read in FASTA sequences and find Ngaps

while (<IN>) {
	chomp;
	if (/^>([^\s]+)/) { 
		$ngaps{$seq_id} = &find_ngaps(\$seq);
		$seq = '';
		$seq_id = $1;
		next;
	}
	$seq .= $_;
}
$ngaps{$seq_id} = &find_ngaps(\$seq);
$seq = '';
$seq_id = '';
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
    my $feature_table  = $doc->createAndAddFeatureTable($seq_stub);
	foreach my $ngap_ref(@{$ngaps{$seq_id}}) {
		my $ngap = $doc->createAndAddFeature(
											$feature_table, 
                                        	$idcreator->new_id(
																db => $options{'project'},
                                             					so_type => 'gap',
																prefix  => $options{'command_id'}
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
    ## check that input and output file parameters were provided	
	unless (defined($options{'input'}) && defined($options{'output'})) {
			$logger->logdie("--input and --output are required parameters");
	}

    ## make sure input file exists
    if (! -e $options{'input'}) { $logger->logdie("input file $options{'input'} does not exist") }
    
    ## make sure output file doesn't exist yet
    if (-e $options{'output'}) { $logger->logdie("can't create $options{'output'} because it already exists") }
    
    return 1;
}
