#!/usr/bin/env perl
=head1 NAME

filter_lgt_best_hit.pl - Wrapper script to find best blast result hit per query ID

=head1 SYNOPSIS

 USAGE: filter_lgt_best_hit.pl
       --input_file=/path/to/some/blast.m8
       --output_file=/path/to/transterm.txt
	   --tmp_dir=/tmp
	   --donor_lineage="Wolbachia"
	   --host_lineage="Drosophila"
	   --tax_id_file=/path/to/tax_ids.txt
	   --nodes_file=/path/to/nodes.txt
	   --names_file\/path/to/names.txt
	 [
	   --filter_lineage='vector'
	   --filter_min_overlap=50
	   --trace_mapping=/path/to/trace_file.txt
	   --host=revan.igs.umaryland.edu
	   --db=gi2taxon
	   --collection=gi2taxonnuc
	   --log=/path/to/file.log
       --debug=3
       --help
     ]

=head1 OPTIONS

B<--input_file,-i>
	A -m8 formatted BLASTN results file

B<--output_dir,-o>
	Path name to output directory.

B<--tmp_dir,-T>
	Directory to store temporary files

B<--trace_mapping,-t>
	Optional. Tab-delimited mapping file that lists the trace ID, clone/mate_pair template ID, and directionality of the read

B<--donor_lineage,-d>
	Name of the donor to find LCA on

B<--host_lineage,-h>
	Name of the host/recipient to find LCA on

B<--filter_lineage,-f>
	Optional. Lineage to filter out as bad (example would be 'vector', which would filter out anything that has a best hit to vector or where a hit overlaps vector by more than the filter_min_overlap).

B<--filter_min_overlap,-m>
	Optional.  Minimum overlap length between a trace/read hit and the filter_lineage. More than this and a clone/matepair will be filtered out. Only applicable if filter_lineage is set.  Default is 50 if unset but applicable

B<--tax_id_file>
	Dump file mapping GI accessions to NCBI taxonomy hits

B<--nodes_file>
	Dump file containing information about NCBI Taxonomy ID nodes and their parent/children relationships.  If not provided, the Entrez database will be used instead to determine taxonomy lineage

B<--names_file>
	Dump file containing information about organism names for given NCBI Taxonomy IDs. If not provided, the Entrez database will be used instead to determine taxonomy lineage

B<--tax_id_file>
	Dump file mapping nucleotide GI accessions to NCBI Taxonomy IDs. If not provided, then the MongoDB database will have to gradually built with an NCBI ESummary taxon lookup

B<--host>
	Optional. MongoDB host server

B<--db>
	Optional. The GI-to-taxon database name (default 'gi2taxon'a)

B<--collection>
	Optional. The table that houses the GI-to-taxon information (default 'gi2taxonnuc')

B<--log,-l>
    Logfile.

B<--debug>
    1,2 or 3. Higher values more verbose.

B<--help>
    Print this message

=head1  DESCRIPTION

 DESCRIPTION

=head1  INPUT

    BlastN results that have been formatted into the m8 format

=head1 OUTPUT

    Three output files containing information about the best hit and LCA for a given query.
	One output will correspond to traces mapping to the host genome
	One output will correspond to traces mapping to the donor genome
	The final output will correspond to the overall mate clones

=head1  CONTACT

    Shaun Adkins
    sadkins@som.umaryland.edu

=cut

use strict;
use warnings;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Pod::Usage;
use LGT::LGTsam2lca;
use LGT::LGTBestBlast;
use GiTaxon;

############# GLOBALS AND CONSTANTS ################
my $debug = 1;
my ($ERROR, $WARN, $DEBUG) = (1,2,3);
my $logfh;
my $no_flatfiles = 0;

my $HOST = 'revan.igs.umaryland.edu:10001';
my $DB = 'gi2taxon';
my $COLL = 'gi2taxonnuc';
my $SAMTOOLS_BIN = '/usr/local/bin/samtools';
####################################################

my %options;

my $results = GetOptions (\%options,
                         "input_file|i=s",
                         "output_dir|o=s",
						 'tmp_dir|T=s',
						 'trace_mapping|t=s',
						 'donor_lineage|d=s',
						 'host_lineage|h=s',
						 'filter_lineage|f=s',
						 'filter_min_overlap|m=s',
						 'tax_id_file=s',
						 'nodes_file=s',
						 'names_file=s',
						 'host|h=s',
						 'db|d=s',
						 'collection|c=s',
                         "log|l=s",
                         "debug=s",
                         "help"
                          );

&check_options(\%options);


$options{host} = $HOST if (! $options{host});
$options{db} = $DB if (! $options{db});
$options{collection} = $COLL if (! $options{collection});

my $gi_tax_obj = GiTaxon->new({
		'nodes' 		=> $options{nodes_file},
		'names' 		=> $options{names_file},
		'gi2tax'		=> $options{tax_id_file},
		'host'			=> $options{host},
		'gi_db'			=> $options{db},
		'gi_coll'		=> $options{collection},
		'taxonomy_dir'	=> $options{tmp_dir},
        'no_flatfiles'  => $no_flatfiles,
		'verbose'		=> 1
	});

my $files = LGT::LGTBestBlast->filterBlast({
		'gitaxon'				=> $gi_tax_obj,
		'output_dir'			=> $options{output_dir},
		'blast'					=> $options{input_file},
		'lineage1'				=> $options{donor_lineage},
		'lineage2'				=> $options{host_lineage},
		'trace_mapping'			=> $options{trace_mapping},
		'filter_lineage'		=> $options{filter_lineage},
		'filter_min_overlap'	=> $options{filter_min_overlap}
	});

exit(0);

sub check_options {
   my $opts = shift;
   if( $opts->{'help'} ) {
       &_pod;
   }

   if( $opts->{'log'} ) {
       open( $logfh, "> $opts->{'log'}") or die("Can't open log file ($!)");
   }

   $debug = $opts->{'debug'} if( $opts->{'debug'} );

   foreach my $req ( qw(input_file output_dir) ) {
       &_log($ERROR, "Option $req is required") unless( $opts->{$req} );
   }

   if (!( -e $opts->{'names_file'} && -e $opts->{'nodes_file'})){
       &_log($DEBUG, "Both --names_file and --nodes_file were not provided or invalid.  Will use Entrez to determine taxonomy lineage instead.");
       $no_flatfiles = 1;
   }

   if (!( -e $opts->{'tax_id_file'})) {
       &_log($DEBUG, "The --tax_id_file option was not passed or invalid.  Will build MongoDB database with NCBI ESummary taxon lookup.");
       $opts->{'tax_id_file'} = '';
   }
}

sub _log {
   my ($level, $msg) = @_;
   if( $level <= $debug ) {
      print STDOUT "$msg\n";
   }
   print $logfh "$msg\n" if( defined( $logfh ) );
   exit(1) if( $level == $ERROR );
}

sub _pod {
    pod2usage( {-exitval => 0, -verbose => 2, -output => \*STDERR} );
}
