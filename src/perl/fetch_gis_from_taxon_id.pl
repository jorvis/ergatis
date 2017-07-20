#!/usr/bin/perl -w

=head1 NAME

fetch_gis_from_taxon_id.pl - Fetch GIs given a taxon id 

=head1 SYNOPSIS

USAGE: fetch_gis_from_taxon_id.pl 
    --database=name of the database to be searched
    --taxon_id=<1234>
    --format=output format
    --output_dir=/path/to/somedir	
=head1 OPTIONS

B<--database,-d>
The name of the database to be searched 

B<--output_dir,-o>
The directory to which the output files will be written.

B<--taxon_id, -t>
The highest-ranking taxon ID to which retrieved GIs would fall under.

B<--format,-f>
The format in which output is required. Examples fasta, gb, asn. If format given is not recognized, eutils will return a GI list.

B<--debug,-d> 
Debug level.  Use a large number to turn on verbose debugging. 

B<--log,-l> 
Log file

B<--help,-h>
This help message

=head1  DESCRIPTION

This script is used to fetch a list of UIs given a taxon ID in a given database.

The script is a modification of fetch_genbank.pl

=head1  CONTACT

Shaun Adkins
sadkins@som.umaryland.edu

=cut

# Define library for the 'get' function used in the next section.
# $utils contains route for the utilities.

use strict;
use lib "/local/projects/ergatis/package-nightly/lib/perl5";
use LWP::Simple;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Pod::Usage;
use POSIX;
use Ergatis::Logger;


my %options = ();
my $results = GetOptions (\%options,
		'output_dir|o=s',
		'database|d=s',
		'taxon_id|t=i',
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

my $id = $options{taxon_id};
chomp($id);

# Remove leading and trailing quotes (since we don't want to convert them later)
$id =~ s/^\"//g;
$id =~ s/\"$//g;

# Substitute spaces for plus signs (to make compatible with GET)
$id=~s/\s+/+/g;
# Substitute quotes and number signs
$id=~s/\"/%22/g;
$id=~s/#/%23/g;

# Let Entrez know this is the organism
my $query = 'txid' . $id . "[ORGN]";

my $ext = defined $options{'format'} ? $options{'format'} : ".gilist";

## open the output file for writing 
my $filename = $options{'output_dir'}."/$id.$ext"; 

open(my $fh, ">$filename") || $logger->logdie("Could not open $filename file for writing");

# ---------------------------------------------------------------------------
# This part of the code is taken from efetch.pl from the NCBI Efetch utilities
# $esearch contains the PATH & parameters for the ESearch call
# $esearch_result containts the result of the ESearch call
# the results are displayed Ánd parsed into variables 
# $Count, $QueryKey, and $WebEnv for later use

my $esearch = "$utils/esearch.fcgi?" .
	"db=$options{'database'}&usehistory=y&term=";

#print "GET $esearch$id\n";

my $esearch_result = get($esearch . $id);

#parse WebEnv, QueryKey and Count (# records retrieved)
my $web = $1 if ($esearch_result =~ /<WebEnv>(\S+)<\/WebEnv>/);
my $key = $1 if ($esearch_result =~ /<QueryKey>(\d+)<\/QueryKey>/);
my $count = $1 if ($esearch_result =~ /<Count>(\d+)<\/Count>/);

die("Had issues getting esearch_result") if (! $count || ! $key || ! $web);


# For larger datasets, EUtils has to download in batches
my $retmax=1000;
for(my $retstart = 0; $retstart < $count; $retstart += $retmax) {
	if( $retstart > 0 ) {
	# sleep in between consecutive requests
	# but not the first time
	#sleep 2;
	}

	my $efetch = "$utils/efetch.fcgi?" .
	"rettype=$options{'format'}&retmode=text&retstart=$retstart&retmax=$retmax&" .
	"db=$options{'database'}&query_key=$key&WebEnv=$web";  
	
	#print "GET $efetch\n";

	my $efetch_result = get($efetch);
	print $fh $efetch_result;
}

close($fh);
# -------------------------------------------------------------------------------

# Subroutine to check the supplied paramaters are correct
sub check_parameters {
	my $options = shift;

	## make sure output_dir, database, query and format were passed
	my @missing;
	foreach my $req( qw(database output_dir query) ) {
	    push(@missing, $req) unless( exists( $options->{$req} ) );
	}
	$logger->logdie("Missing mandatory parameters: ".join(", ", @missing)) unless( @missing == 0 );

	$options->{'format'} = 'gi' if ! $options->{'format'};

	## make sure the output_dir exists, otherwise make it.
	system("mkdir -p $options{'output_dir'}") unless( -d $options{'output_dir'} );
}

