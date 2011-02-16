#!/usr/bin/perl 

eval 'exec /usr/bin/perl  -S $0 ${1+"$@"}'
if 0; # not running under some shell

BEGIN{foreach (@INC) {s/\/usr\/local\/packages/\/local\/platform/}};
use lib (@INC,$ENV{"PERL_MOD_DIR"});
no lib "$ENV{PERL_MOD_DIR}/i686-linux";
no lib ".";

=head1 NAME

fetch_genbank.pl - Fetches data from GenBank for the provided accession numbers list 

=head1 SYNOPSIS

USAGE: fetch_genbank.pl 
--database=name of the database to be searched
--query=comma-separated accession number list
--format=output format
--output_dir=/path/to/somedir	
=head1 OPTIONS

B<--database,-d>
The name of the database to be searched 

B<--output_dir,-o>
The directory to which the output file will be written.

B<--query,-q>
The comma-separated list of accession numbers from the genBank to be fetched.

B<--format,-f>
The format in which output is required like FASTA

B<--output_file,-s>
The name of the output file containing the searched, reference sequence

B<--debug,-d> 
Debug level.  Use a large number to turn on verbose debugging. 

B<--log,-l> 
Log file

B<--help,-h>
This help message

=head1  DESCRIPTION

This script is used to fetch a text query given a UID in a given database. Its mainly used for downloading reference genomes in FASTA format from genBank
given genBank accesion ids. 

=head1  INPUT

=head1  CONTACT

Sonia Agrawal
sagrawal@som.umaryland.edu

=cut

# Define library for the 'get' function used in the next section.
# $utils contains route for the utilities.

use strict;
use LWP::Simple;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Pod::Usage;
use POSIX;
BEGIN {
	use Ergatis::Logger;
}


my %options = ();
my $results = GetOptions (\%options,
		'output_dir|o=s',
		'database|d=s',
		'query|q=s',
		'format|f=s',
		'log|l=s',
		'debug|b=s',
		'help|h') || pod2usage();

my $utils = "http://www.ncbi.nlm.nih.gov/entrez/eutils";

my $logfile = $options{'log'} || Ergatis::Logger::get_default_logfilename();
my $logger = new Ergatis::Logger('LOG_FILE'=>$logfile,
		'LOG_LEVEL'=>$options{'debug'});
$logger = $logger->get_logger();

## display documentation
if( $options{'help'} ){
	pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDERR} );
}

## make sure everything passed was correct
&check_parameters(\%options);
my $report = $options{format};
my $db = $options{database};
my @ids = split(/,/,$options{query}); 
foreach my $id (@ids) {
	chomp($id);
## open the output file for writing 
	my $fh;
	my $filename = $options{output_dir}."/reference_".$id.".".$options{format}; 

	open($fh, ">$filename") || $logger->logdie("Could not open $filename file for writing");

# ---------------------------------------------------------------------------
# This part of the code is taken from efetch.pl from the NCBI Efetch utilities
# $esearch contains the PATH & parameters for the ESearch call
# $esearch_result containts the result of the ESearch call
# the results are displayed Ánd parsed into variables 
# $Count, $QueryKey, and $WebEnv for later use

	my $esearch = "$utils/esearch.fcgi?" .
		"db=$db&retmax=1&usehistory=y&term=";

	my $esearch_result = get($esearch . $id);

#print "\nESEARCH RESULT: $esearch_result\n";

	$esearch_result =~ 
		m|<Count>(\d+)</Count>.*<QueryKey>(\d+)</QueryKey>.*<WebEnv>(\S+)</WebEnv>|s;

	my $Count    = $1;
	my $QueryKey = $2;
	my $WebEnv   = $3;

#print "Count = $Count; QueryKey = $QueryKey; WebEnv = $WebEnv\n";

	my $retstart;
	my $retmax=3;
	for($retstart = 0; $retstart < $Count; $retstart += $retmax) {
		my $efetch = "$utils/efetch.fcgi?" .
			"rettype=$report&retmode=text&retstart=$retstart&retmax=$retmax&" .
			"db=$db&query_key=$QueryKey&WebEnv=$WebEnv";

#print "\nEF_QUERY=$efetch\n";     

		my $efetch_result = get($efetch);

		sleep 2;

		print $fh $efetch_result;
	}

}
# -------------------------------------------------------------------------------

# Subroutine to check the supplied paramaters are correct
sub check_parameters {
	my $options = shift;

## make sure output_dir, database, query and format were passed
	unless ( $options{database} && $options{output_dir} && $options{format} && $options{query}) {
		$logger->logdie("All the manadatory parameters should be passed");
	}

## make sure the output_dir exists
	if (! -e "$options{output_dir}") {
		$logger->logdie("The output directory passed could not be read or does not exist");
	}
}

