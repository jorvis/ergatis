#!/usr/local/packages/perl-5.8.5/bin/perl
use lib (@INC,$ENV{"PERL_MOD_DIR"});
no lib "$ENV{PERL_MOD_DIR}/i686-linux";
no lib ".";

=head1  NAME 

tmhmm2bsml.pl - convert TMHMM output to BSML

=head1 SYNOPSIS

USAGE: tmhmm2bsml.pl --input=/path/to/tmhmm_file --output=/path/to/output.bsml

=head1 OPTIONS

B<--input,-i> 
    Input file file from a tmhmm run [long format only!].

B<--debug,-d> 
    Debug level.  Use a large number to turn on verbose debugging. 

B<--log,-l> 
    Log file

B<--output,-o> 
    Output BSML file (will be created, must not exist)

B<--help,-h> 
    This help message

=head1   DESCRIPTION

This script is used to convert the output from TMHMM into BSML.

=head1 INPUT

TMHMM can be run using multiple input sequences simultaneously, and this
script supports parsing single or multiple input result sets.
 
You define the input file using the --input option.  This file does not need any
special file extension.

=head1 OUTPUT

After parsing the input file, a file specified by the --output option is created.  This script
will fail if it already exists.

=head1 CONTACT

Brett Whitty
bwhitty@tigr.org

=cut

use strict;
use Log::Log4perl qw(get_logger);
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use BSML::BsmlBuilder;
use BSML::BsmlParserTwig;
use BSML::BsmlRepository;
use Papyrus::TempIdCreator;
use Pod::Usage;
use Workflow::Logger;

my %options = ();
my $results = GetOptions (\%options, 
			  'input|i=s',
              'output|o=s',
              'debug|d=s',
              'command_id=s',       ## passed by workflow
              'logconf=s',          ## passed by workflow (not used)
              'project|p=s',
              'log|l=s',
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

## we want to create ids unique to this document, which will be replaced later.  they must
##  contain the prefix that will be used to look up a real id, such as ath1.gen.15
my $next_id = 1;

## we want a new doc
my $doc = new BSML::BsmlBuilder();

## we're going to generate ids
my $idcreator = new Papyrus::TempIdCreator();

## recognized output result data term mappings
my %result_term = (
					'Length' 					=> 'length',
					'Number of predicted TMHs'	=> 'tmh_count',
					'Exp number of AAs in TMHs' => 'exp_aa_in_tmh',
					'Exp number, first 60 AAs'  => 'exp_first_60',
					'Total prob of N-in'		=> 'prob_n_in',
				   );
## recognized results qualifier term mappings
my %qualifier = (
					'POSSIBLE N-term signal sequence' => 'n_term_signal',
				);

## open the input file for parsing
open (my $ifh, $options{'input'}) || $logger->logdie("can't open input file for reading");

my @sequence_ids;
my %result_ref_hash;
#my $temp;
my %property_hash;
my %segment_hash;
my %qualifier_hash;
my $skip_flag = 0;
#my $line;
while (my $line = <$ifh>) {
    chomp $line;
	## skip lines until we hit #'s 
	if ($skip_flag && !($line =~ /^#/)) {
			next;
	}
	if ($line =~ /^# ([^ ]+) ([^:]+):\s+([^ ]+)$/) {
		$skip_flag = 0;
		my $seq_id = $1;
		my $term = $2;
		my $value = $3;
		if (defined($result_term{$term})) {
			if (!defined($property_hash{$seq_id})) {
				$property_hash{$seq_id} = [];
			}
			push(@{$property_hash{$seq_id}}, [$seq_id, $result_term{$term}, $value]);
			## We may want to use these in the output later
			#print STDERR $result_term{$term};
		} else {
			die "unrecognized result term";
		}
	} elsif ($line =~ /^# ([^ ]+) (.*)$/) {
		$skip_flag = 0;
		my $seq_id = $1;
		my $term = $2;
		if (defined($qualifier{$term})) {
			## We may want to use these in the output later
			#print STDERR $qualifier{$term}."\n";
			if (!defined($qualifier_hash{$seq_id})) {
				$qualifier_hash{$seq_id} = {};
			}
			#push(@{$qualifier_hash{$seq_id}}, [$seq_id, $result_term{$term}, $value]);
			$qualifier_hash{$seq_id}->{$qualifier{$term}} = 1;
		} else {
			die "unrecognized qualifier term";
		}
	}
	if ($line =~ /^([^\s]+)\s+([^\s]+)\s+([^\s]+)\s+(\d+)\s+(\d+)/) {
		my $id = $1;
		my $model = $2;
		my $state = $3;
		my $start = $4;
		my $end = $5;
		if (!defined($segment_hash{$id})) {
			$segment_hash{$id} = [];
		}
		push(@{$segment_hash{$id}}, [$id, $model, $state, $start, $end]);
	}
}

#print $segment_hash{$id}->[0]->[0];

#die();

foreach my $s(keys(%segment_hash)) {
   	my $seq = $doc->createAndAddSequence($s, undef, '', 'aa', 'polypeptide');
   	$seq->addBsmlLink('analysis', '#tmhmm_analysis', 'input_of');
	
	foreach my $prop_ref(@{$property_hash{$s}}) {
		$doc->createAndAddBsmlAttribute($seq, $prop_ref->[1], $prop_ref->[2]);
	}
	
	my $feature_table;
	$feature_table  = $doc->createAndAddFeatureTable($seq);
								   
	if ($qualifier_hash{$s}->{'n_term_signal'}) {							   
		my $new_id = $idcreator->new_id( db      => $options{project},
			                             so_type => 'signal_peptide',
										 prefix  => $options{command_id},
									   );
		my $signalp = $doc->createAndAddFeature($feature_table, 
			                                    $new_id, 
												'', 
												'signal_peptide'
											   );
										   
	 	$signalp->addBsmlLink('analysis', '#tmhmm_analysis', 'computed_by');
	}
	
	my $new_id = $idcreator->new_id( db      => $options{project},
       	                             so_type => 'located_sequence_feature',
									 prefix  => $options{command_id},
               	                   );

	foreach my $tm_ref(@{$segment_hash{$s}}) {
		my $tm = $doc->createAndAddFeature(
											$feature_table, 
                                        	$idcreator->new_id(
																db => $options{'project'},
                                             					so_type => 'located_sequence_feature',
																prefix  => $options{'command_id'}
															  ),
                                          	$tm_ref->[2], ## state
											'located_sequence_feature'
										);
                                          
	    $tm->addBsmlLink('analysis', '#tmhmm_analysis', 'computed_by');
    	$tm->addBsmlIntervalLoc(
									$tm_ref->[3] - 1, ## start
								   	$tm_ref->[4], ## end
								   	0
								 );
	}

}
    
## add the analysis element
	my $analysis = $doc->createAndAddAnalysis(
                            id => 'tmhmm_analysis',
                            sourcename => $options{'output'},
                          );

## now write the doc
$doc->write($options{'output'});

exit();

sub check_parameters {
    
    ## make sure input file exists
    if (! -e $options{'input'}) { $logger->logdie("input file $options{'input'} does not exist") }
    
    ## make sure output file doesn't exist yet
    if (-e $options{'output'}) { $logger->logdie("can't create $options{'output'} because it already exists") }
    
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

