#!/usr/local/bin/perl

use lib (@INC,$ENV{"PERL_MOD_DIR"});
no lib "$ENV{PERL_MOD_DIR}/i686-linux";
no lib ".";

=head1  NAME 

pepstats2bsml.pl - convert EMBOSS pepstats output to BSML

=head1 SYNOPSIS

USAGE: pepstats2bsml.pl --input=/path/to/pepstats_out.raw --output=/path/to/output.bsml

=head1 OPTIONS

B<--input,-i> 
    Input file file from a pepstats run.

B<--debug,-d> 
    Debug level.  Use a large number to turn on verbose debugging. 

B<--log,-l> 
    Log file

B<--output,-o> 
    Output BSML file (will be created, must not exist)

B<--help,-h> 
    This help message

=head1   DESCRIPTION

This script is used to convert the output from pepstats into BSML.

=head1 INPUT

Pepstats can be run using multiple input sequences simultaneously, but this
script supports parsing output from single sequence input files only. If it
is run on output from multiple sequences, it will only parse the first entry.
 
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

## map output to ontology terms in output.obo
my $onto_terms = {
'Molecular weight' 									=> 'mol_wt',
'Residues' 											=> '',
'Average Residue Weight' 							=> 'avg_residue_wt',
'Charge' 											=> 'charge',
'Isoelectric Point' 								=> 'isoelectric_pt',
'A280 Molar Extinction Coefficient' 				=> 'extinction_coefficient_mol',
'A280 Extinction Coefficient 1mg/ml' 				=> 'extinction_coefficient_mg_ml',
'Improbability of expression in inclusion bodies' 	=> 'improb_expr_in_inclusion_bodies',
'Probability of expression in inclusion bodies' 	=> 'prob_expr_in_inclusion_bodies',
		    	 };

my $doc = new BSML::BsmlBuilder();
$doc->makeCurrentDocument();

my $seq_id;
my $atts;
my $length;

open (IN, $options{'input'}) || $logger->logdie("couldn't open input file for reading");
while (<IN>) {
	chomp;
	s/\t+/\t/g;
	s/\s+=\s+/=/g;
	s/\s+$//;
	if (/^\s*$/) {
		next;
	}
	if (/PEPSTATS of ([^\s]+)/) {
		$seq_id = $1;
		next;
	}
	if (/^Residue/) {
		last;
	}
	my @data = split("\t");
	foreach my $datum(@data) {
		my ($att, $val) = split("=",$datum);
		if ($att eq 'Residues') {
			$length = $val;
		}
		if (!defined($onto_terms->{$att})) { 
			$logger->logdie("encountered an unexpected data type '$att'");
		} elsif ($onto_terms->{$att} ne '') {
			$atts->{$onto_terms->{$att}} = $val;
		}
	}
}
my $seq = $doc->createAndAddSequence(
										$seq_id, 
										$seq_id, 
										$length, 
										'aa', 
										'polypeptide'
									);
foreach my $att(keys(%{$atts})) {
	$seq->addBsmlAttr($att, $atts->{$att});
}
$doc->write($options{'output'});
