#!/usr/bin/env perl 

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
The directory to which the output files will be written.

B<--query,-q>
The comma-separated list of accession numbers from the genBank to be fetched.

B<--format,-f>
The format in which output is required. Examples fasta, gb, asn. If format given is not recognized, eutils will return asn.

B<--debug,-d> 
Debug level.  Use a large number to turn on verbose debugging. 

B<--log,-l> 
Log file

B<--help,-h>
This help message

=head1  DESCRIPTION

This script is used to fetch a text query given a UID in a given database. Its mainly used for downloading reference genomes in FASTA format from genBank
given genBank accesion ids. 

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
use Ergatis::Logger;


my %options = ();
my $results = GetOptions (\%options,
		'output_dir|o=s',
		'database|d=s',
		'query|q:s',
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

foreach my $id ( split(/,/, $options{query}) ) {
	chomp($id);
	next if ($id =~ /\/|\\/g);

	## open the output file for writing 
	my $filename = $options{'output_dir'}."/".$id.".".$options{'format'}; 

	open(my $fh, ">$filename") || $logger->logdie("Could not open $filename file for writing");

	# ---------------------------------------------------------------------------
	# This part of the code is taken from efetch.pl from the NCBI Efetch utilities
	# $esearch contains the PATH & parameters for the ESearch call
	# $esearch_result containts the result of the ESearch call
	# the results are displayed Ánd parsed into variables 
	# $Count, $QueryKey, and $WebEnv for later use

	my $esearch = "$utils/esearch.fcgi?" .
		"db=$options{'database'}&retmax=1&usehistory=y&term=";

	my $esearch_result = get($esearch . $id);

	my ($Count, $Id);
	my $retmax=3;
	if( $esearch_result =~ m|<Count>(\d+)</Count>.*<Id>(\d+)</Id>|s ) {
	    $Count = $1;
	    $Id = $2;
	} else {
	    die("Had issues getting esearch_result");
	}

	for(my $retstart = 0; $retstart < $Count; $retstart += $retmax) {
	    if( $retstart > 0 ) {
		# sleep in between consecutive requests
		# but not the first time
		sleep 2;
	    }

	    my $efetch = "$utils/efetch.fcgi?" .
		"rettype=$options{'format'}&retmode=text&retstart=$retstart&retmax=$retmax&" .
		"db=$options{'database'}&id=$Id";  

	    my $efetch_result = get($efetch);
	    print $fh $efetch_result;
	}

	close($fh);

}
# -------------------------------------------------------------------------------

# Subroutine to check the supplied paramaters are correct
sub check_parameters {
	my $options = shift;

	## make sure output_dir, database, query and format were passed
	my @missing;
	foreach my $req( qw(database output_dir format query) ) {
	    push(@missing, $req) unless( exists( $options->{$req} ) );
	}
	$logger->logdie("Missing mandatory parameters: ".join(", ", @missing)) unless( @missing == 0 );

	## make sure the output_dir exists, otherwise make it.
	system("mkdir -p $options{'output_dir'}") unless( -d $options{'output_dir'} );
}

