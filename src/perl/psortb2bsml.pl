#!/usr/bin/perl
=head1  NAME 

psortb2bsml.pl - convert psortb output to BSML

=head1 SYNOPSIS

USAGE: psortb2bsml.pl --input=/path/to/psortb_out.raw --output=/path/to/output.bsml

=head1 OPTIONS

B<--input,-i> 
    Input file file from a psortb run.

B<--debug,-d> 
    Debug level.  Use a large number to turn on verbose debugging. 

B<--log,-l> 
    Log file

B<--output,-o> 
    Output BSML file (will be created, must not exist)

B<--help,-h> 
    This help message

=head1   DESCRIPTION

This script is used to convert the output from psortb into BSML.

=head1 INPUT

Pepstats can be run using multiple input sequences simultaneously, and this
script will support this.
 
You define the input file using the --input option.  This file does not need any
special file extension.

=head1 OUTPUT

After parsing the input file, a file specified by the --output option is 
created.  This script will fail if it already exists.

=head1 CONTACT

Brett Whitty
bwhitty@tigr.org

=cut

use Log::Log4perl qw(get_logger);
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use Pod::Usage;
use Workflow::Logger;
use BSML::BsmlBuilder;
use warnings;
use strict;

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

if (!$options{'input'}) {
	pod2usage("you must specify an input file with --input");
}
if (!$options{'output'}) {
	pod2usage("you must specify an output file with --output");
}
if (-e $options{'output'}) {
	$logger->logdie("output file '$options{output}' already exists");
}
if (!-e $options{'input'}) {
	$logger->logdie("input file '$options{input}' does not exist");
}

## ontology terms in output.obo for the output columns
my @columns = (
	'cmsvm_loc',
	'cmsvm_desc',
	'cytosvm_loc',
	'cytosvm_desc',
	'ecsvm_loc',
	'ecsvm_desc',
	'hmmtop_loc',
	'hmmtop_desc',
	'motif_loc',
	'motif_desc',
	'ompmotif_loc',
	'ompmotif_desc',
	'omsvm_loc',
	'omsvm_desc',
	'ppsvm_loc',
	'ppsvm_desc',
	'profile_loc',
	'profile_desc',
	'scl-blast_loc',
	'scl-blast_desc',
	'scl-blaste_loc',
	'scl-blaste_desc',
	'signal_loc',
	'signal_desc',
	'cytoplasmic_score',
	'cytoplasmic_membrane_score',
	'periplasmic_score',
	'outer_membrane_score',
	'extracellular_score',
	'final_loc',
	'final_score',
);

my $doc = new BSML::BsmlBuilder();
$doc->makeCurrentDocument();

my $seq_id;
my $atts;
my $length;

open (IN, $options{'input'}) || $logger->logdie("couldn't open input file for reading");
my $header = <IN>; ## skip the header line
if (!($header =~ /^SeqID/)) {
	$logger->logdie("invalid pSortB file or file format has changed");
}

my $psortb;

while (<IN>) {
	chomp;
	if (/^\s*$/) {
		next;
	}
	my @cols = split("\t");
	if (scalar(@cols) != 32) {
		$logger->logdie("expected 32 cols in output but found ".scalar(@cols));
	}
	my $seq_id = shift @cols;
	$seq_id =~ s/^([^\s]+).*/$1/;
	$psortb->{$seq_id} = \@cols;
}

foreach my $seq_id(keys(%{$psortb})) {
	my $seq = $doc->createAndAddSequence(
										$seq_id, 
										$seq_id, 
										'', 
										'aa', 
										'polypeptide'
									);
	for (my $i=0; $i<31; $i++) { 
		$seq->addBsmlAttr($columns[$i], $psortb->{$seq_id}->[$i]);
	}
	$doc->createAndAddLink(
						$seq,
                        'analysis',
                        '#psortb_analysis',
                        'computed_by',
                      );



}
$doc->write($options{'output'});
