#!/usr/bin/env perl
=head1 NAME

blast2lca.pl - Wrapper script to generate Lowest Common Ancestors from a BAM/SAM file

=head1 SYNOPSIS

 USAGE: blast2lca.pl
       --input_file=/path/to/some/input.bam
       --output_dir=/path/to/transterm
	   --tmp_dir=/tmp
	   --tax_id_file=/path/to/tax_ids.txt
	   --nodes_file=/path/to/nodes.txt
	   --names_file\/path/to/names.txt
	 [
	   --host=revan.igs.umaryland.edu
	   --db=gi2taxon
	   --collection=gi2taxonnuc
	   --log=/path/to/file.log
       --debug=3
       --help
     ]

=head1 OPTIONS

B<--input_file,-i>
	An -m8 formatted blast results file

B<--output_directory,-o>
	Path name to LCA output directory.

B<--evalue_cutoff, -e>
	Max e-value allowed for a hit. Default is '1' if not specified

B<--best_hits_only, -b>
	Parse the Blast results for only the best hit per query, as created by lgt_best_blast_hit.pl
    Default is 1

B<--tmp_dir,-t>
	Directory to store temporary files

B<--tax_id_file>
	Dump file mapping GI accessions to NCBI taxonomy IDs

B<--nodes_file>
    Dump file containing information about NCBI Taxonomy ID nodes and their parent/children relationships.  If not provided, the Entrez database will be used instead to determine taxonomy lineage

B<--names_file>
    Dump file containing information about organism names for given NCBI Taxonomy IDs. If not provided, the Entrez database will be used instead to determine taxonomy lineage

B<--tax_id_file>
	Dump file mapping nucleotide GI accessions to NCBI Taxonomy IDs. If not provided, then the MongoDB database will have to gradually built with an NCBI ESummary taxon lookup

B<--host>
	MongoDB host server

B<--db>
	The GI-to-taxon database name (default 'gi2taxon')

B<--collection>
	The table that houses the GI-to-taxon information (default 'gi2taxonnuc')

B<--log,-l>
    Logfile.

B<--debug,-D>
    1,2 or 3. Higher values more verbose.

B<--help,-H>
    Print this message

=head1  DESCRIPTION

 DESCRIPTION

=head1  INPUT

    Describe the input

=head1 OUTPUT

    Two LCA output files.

=head1  CONTACT

    Shaun Adkins
    sadkins@som.umaryland.edu

=cut

use strict;
use warnings;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev pass_through);
use Pod::Usage;
use LGT::LGTSeek;
use GiTaxon;

############# GLOBALS AND CONSTANTS ################
my $debug = 1;
my ($ERROR, $WARN, $DEBUG) = (1,2,3);
my $logfh;
my $no_flatfiles = 0;

my $HOST = 'revan.igs.umaryland.edu:10001';
my $DB = 'gi2taxon';
my $COLL = 'gi2taxonnuc';
my $EVAL_CUTOFF = 1;
my $BEST_HITS_FLAG = 1;
####################################################

my %options;

my $results = GetOptions (\%options,
                         "input_file|i=s",
                         "output_dir|o=s",
						 'tmp_dir|T=s',
						 'tax_id_file=s',
						 'nodes_file=s',
						 'names_file=s',
						 'host|h=s',
						 'db|d=s',
						 'collection|c=s',
                         'evalue_cutoff|e=s',
						 'best_hits_only|b=i',
						 "log|l=s",
                         "debug|D=i",
                         "help|H"
                          );

&check_options(\%options);

$options{host} = $HOST if (! $options{host});
$options{db} = $DB if (! $options{db});
$options{collection} = $COLL if (! $options{collection});
$options{evalue_cutoff} = $EVAL_CUTOFF if (! $options{evalue_cutoff});
$options{best_hits_flag} = $BEST_HITS_FLAG if (! $options{best_hits_flag});

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
my $lgtseek_obj = LGT::LGTSeek->new({
		'verbose'		=> 1,
		'output_dir'		=> $options{output_dir},
	});

my $files = $lgtseek_obj->blast2lca({
		'gi2tax'			=> $gi_tax_obj,
		'blast'				=> $options{input_file},
		'evalue_cutoff'		=> $options{evalue_cutoff},
		'best_hits_only'	=> $options{best_hits_only},
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
