#!/local/perl/bin/perl

eval 'exec /local/perl/bin/perl  -S $0 ${1+"$@"}'
    if 0; # not running under some shell

=head1 NAME

db2bsml - Dumps a database to a BSML document(s)

=head1 SYNOPSIS

USAGE:  db2bsml -u username -p password -d databases [-a assembly_name] [-b organism] [-o output_dir] [-t dtd] [-s schema] [-v]

=head1 OPTIONS

=over 8

=item B<--username,-u>
    
    Database username

=item B<--password,-p>
    
    Database password

=item B<--database,-d>
    
    Database name

=item B<--output_dir,-o> Directory to save output. Default directory
    is pulled from BSML_REPOSITORY environment variable

=item B<--asmbl_ids,-a>

    Optional. User can specify a single assembly name/id or a comma
    separated list.  If this option is ommited all asmbl_ids will be
    dumped.

=item B<--organism,-b>

    Optional. User can specify a organism name/id.  This option
    overrides the --asmbl_id option if they are both specified

=item B<--dtd,-t>

    Optional. DTD to validate against.  If blank validate to inline DTD.

=item B<--schema,-s>

    Optional. Schema to validate against. If blank validate to inline schema

=item B<--log,-l>
    
    Optional.  Write debug information to a log file.

=item B<--help,-h>

    Print this help


=back

=head1 DESCRIPTION

    db2bsml dumps database content to BSML.  One document is created
    per assembly in the output directory.

=cut

use lib "shared";
use strict;
use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use Pod::Usage;
use BSML::BsmlBuilder;

#Process command line options
my ($genus, $species, $strain, $organism, $source_database, $fasta_file, $output_file, $help);
my $results = GetOptions (
			  'help|h=s'        => \$help,
			  'genus|g=s'       => \$genus,
			  'species|s=s'     => \$species,
			  'strain|S=s'      => \$strain,
			  'organism|o=s'    => \$organism,
			  'output_file|o=s'  => \$output_file,
			  'fasta_file|f=s'  => \$fasta_file,
			  'source_database|d=s' => \$source_database
			  );

# pod2usage({-exitval => 1, -verbose => 2, -output => \*STDOUT}) if($help || !$genus || !$species || !$source_database || !$fasta_file || !$output_file);

my $doc = new BSML::BsmlBuilder(); 

$organism = $doc->createAndAddOrganism( 
					'genome'  => $doc->createAndAddGenome(),
					'genus'   => $genus,
					'species' => $species
					);

$strain = $doc->createAndAddStrain( 
				    'organism'        => $organism,
				    'name'            => $strain, 
				    'database'        => $source_database,
				    'source_database' => $source_database
				    );				   

my ($uid, $seq, $output_dir);
open(F_IN, $fasta_file) || die ("no can find: $fasta_file");
my $first_time = 1;
while(<F_IN>) {
    s/\n//;
    if(/^>/) {
	s/>//;
	s/\s.*//;
	$uid = $_;

	$uid = "_" . $uid if ($uid =~ /^[0-9]/);

	if ($first_time != 1) {
	    add_stuff($doc, $uid, $seq, $fasta_file);
	}

	$first_time = 0;
	my $seq;
    }
    else {
	$seq .= $_;
    }
}
add_stuff($doc, $uid, $seq, $fasta_file);

$doc->write("$output_file");

if(! -e "$output_file"){
    die ("File not created $output_file");
}

sub add_stuff {
    my($doc, $uid, $seq, $fasta_file) = @_;
    my($seq_length);

    $seq_length = length($seq);

    my $asmseq = $doc->createAndAddExtendedSequenceN( 'id' => $uid, 
						      'title' => '', 
						      'length' => $seq_length, 
						      'molecule' => 'aa', 
						      'locus' => '', 
						      'dbsource' => '', 
						      'icAcckey' => '', 
						      'strand' => '');

    $doc->createAndAddSeqDataImport($asmseq, 
				    'fasta', 
				    $fasta_file, 
				    $uid);
}

