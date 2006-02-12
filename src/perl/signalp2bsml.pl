#!/usr/local/packages/perl-5.8.5/bin/perl
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

use strict;
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

## open the input file for parsing
open (my $ifh, $options{'input'}) || $logger->logdie("can't open input file for reading");

my @sequence_ids;
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
	my $sequence_id;
	if ($next_seq_flag) {
		$next_seq_flag = 0;
		$sequence_id = $_;
	} else {
		$sequence_id = <$ifh>;
	}
	chomp $sequence_id;
	$sequence_id =~ s/^>//;
	$logger->debug("sequence id is '$sequence_id'") if($logger->is_debug());  
	push(@sequence_ids, $sequence_id);
		
	while (my $result_line = <$ifh>) {
		chomp $result_line;
		if ($result_line =~ /^SignalP-NN result:/) {
			$logger->debug("parsing signalp nn results") if($logger->is_debug());  
			$temp = <$ifh>;
			chomp $temp;
			$temp =~ /length = (\d+)/ || $logger->error("Couldn't parse sequence length in NN result.") if ($logger->is_error);
#			$nn_results{'length'} = $1;
			$temp = <$ifh>; # we will discard the header line
			$temp = <$ifh>;
			chomp $temp;
			$temp =~ /^  max\. C\s+([^\s]+)\s+([^\s]+)\s+([^\s]+)\s+([^\s]+)/ || $logger->logdie("Couldn't parse max. C line in NN result");
			#$temp =~ /^  max\. C[ ]+([^ ]+)[ ]+([^ ]+)[ ]+([^ ]+)[ ]+([^ ]+)/ || $logger->logdie("Couldn't parse max. C line in NN result");
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
				$result_ref_hash{$sequence_id, 'nn'} = \%nn_results;
				$nn_results{'cleavage_site_position'} = $1."-".$2;
			} else {
				$result_ref_hash{$sequence_id, 'nn'} = \%nn_results;
			}

		} elsif ($result_line =~ /^SignalP-HMM result:/) {
			$logger->debug("parsing signalp hmm results") if($logger->is_debug());  
			$temp = <$ifh>; # we will discard the header line
			$temp = <$ifh>;
			chomp $temp;
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
			#print $temp."\n";
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

    foreach my $s(@sequence_ids) {
    	my $seq = $doc->createAndAddSequence($s, undef, '', 'aa', 'polypeptide');
	foreach my $a(('nn','hmm')) {
	   	$seq->addBsmlLink('analysis', '#signalp_'.$a.'_analysis', 'input_of');
	}
	
	my $ft;
	foreach my $a(('nn','hmm')) {
		if (
			($a eq 'hmm' && $result_ref_hash{$s,$a}->{'prediction'} eq 'Signal peptide' )
			|| 
			($a eq 'nn' && $result_ref_hash{$s,$a}->{'d_signal_peptide'} eq 'YES') ) {
		if (!$ft){
        		$ft  = $doc->createAndAddFeatureTable($seq);
		}
		my $new_id = $idcreator->new_id( db      => $options{project},
                                                 so_type => 'signal_peptide',                                           	 prefix  => $options{command_id},
                                                );
        		my $signalp = $doc->createAndAddFeature($ft, $new_id, '', 'signal_peptide');
	   	   	$signalp->addBsmlLink('analysis', '#signalp_'.$a.'_analysis', 'computed_by');

			foreach my $k(keys %{$result_ref_hash{$s,$a}}) {
				$doc->createAndAddBsmlAttribute($signalp, $k, $result_ref_hash{$s,$a}->{$k} );
			}
		}
	}
    }
    
## add the analysis element
foreach my $a(('nn','hmm')) {
	my $analysis = $doc->createAndAddAnalysis(
                            id => 'signalp_'.$a.'_analysis',
                            sourcename => $options{'output'},
                          );
}

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

