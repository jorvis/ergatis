#! /usr/local/bin/perl

# this is temporary for development
# use lib "/export/CVS/bsml/src/";

=head1  NAME 

Tiling2Bsml.pl  - convert show_tiling output into BSML sequence mapping documents

=head1 SYNOPSIS

USAGE:  Tiling2Bsml.pl -t tilings_file -b path_to_bsml_repository -o output_file_bsml

=head1 OPTIONS

=over 4

=item *

B<--tilingPath,-t> [REQUIRED] file containing tiling path data from show_tiling

=item * 

B<--bsml_repository,-b> [REQUIRED] filepath to the bsml repository

=item *

B<--output,-o> [REQUIRED] output BSML file containing coverage information

=item * 

B<--help,-h> This help message

=back

=head1   DESCRIPTION

Tiling2Bsml.pl is designed to convert tiling path data from show_tiling 
into BSML sequence mapping documents.


The tiling path input should be tab delimited with the following fields. The 
first line encode the reference assembly id, starting with a '>'

# 0 - start in ref
# 1 - end in ref
# 2 - distance to next contig
# 3 - length of this contig
# 4 - alignment coverage
# 5 - identity
# 6 - orientation
# 7 - id

=cut

use strict;
use warnings;

use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use Pod::Usage;

BEGIN {
    require '/usr/local/devel/ANNOTATION/ard/chado-v1r5b1/lib/site_perl/5.8.5/BSML/BsmlBuilder.pm';
    import BSML::BsmlBuilder;
}

my %options = ();
my $results = GetOptions( \%options, 'tilingPath|t=s', 'bsml_repository_path|b=s', 'outFile|o=s', 'help|h' ) || pod2usage();

# display documentation
if( $options{'help'} ){
    pod2usage( {-exitval=>0, -verbose => 2, -output => \*STDOUT} );
}

if( !( $options{'tilingPath' } ) )
{
    pod2usage( {-exitval=>1, -verbose => 2, -output => \*STDOUT} );
}

if( !( $options{'bsml_repository_path' } ) )
{
    pod2usage( {-exitval=>1, -verbose => 2, -output => \*STDOUT} );
}

if( !( $options{'outFile' } ) )
{
    pod2usage( {-exitval=>1, -verbose => 2, -output => \*STDOUT} );
}

# The output of show_tiling is tab deliminted with the following fields
# 0 - start in ref
# 1 - end in ref
# 2 - distance to next contig
# 3 - length of this contig
# 4 - alignment coverage
# 5 - identity
# 6 - orientation
# 7 - id

my $repositoryPath = $options{'bsml_repository_path'};

open( TILINGS, $options{'tilingPath'} ) or die "Unable to open $options{'tilingPath'}\n";

my $bsmlDoc = new BSML::BsmlBuilder;

# add the reference sequence to the Tiling Document

my $reference_assmbl_id = '';

my $line = <TILINGS>;

if( !($line) )
{
    die "Tiling path contains no reference sequence\n";
}
else
{
    # Grab the refseq id from the first line
    if ($line =~ />([\S]*)[\s]([\d]*)/)
    {
	$reference_assmbl_id = $1;

	# add the reference sequence to the bsml tiling
	my $newSeq = $bsmlDoc->createAndAddExtendedSequenceN( id => $1,
							      title => $1,
							      molecule => 'dna',
							      length => $2 );

	$newSeq->addattr( 'class', 'assembly' );

	# add a file link to the assembly bsml doc in the repository

	$bsmlDoc->createAndAddSeqDataImportN( seq => $newSeq,
				      format => "BSML",
                                      source => $repositoryPath."/".$reference_assmbl_id.".bsml",
                                      id => "_$reference_assmbl_id");
	
    }
    else
    {
	die "Unable to determine reference sequence\n";
    }
}

my $rank = 0;

# add each tiled sequence to the Tiling Document

while( my $line = <TILINGS> )
{
    my @tile = split( "\t", $line );
    chomp( $tile[7] );

    my $assmbl_id = $tile[7];
    my $orientation = $tile[6];
    my $ascending = 0;

    if( $orientation eq '+' )
    {
	$ascending = 1;
    }
    else
    {
	$ascending = 0;
    }

    my $newSeq = $bsmlDoc->createAndAddExtendedSequenceN( id => $assmbl_id,
							  title => $assmbl_id,
							  molecule => 'dna',
							  length => $tile[3] );

    $newSeq->addattr( 'class', 'assembly' );

    $bsmlDoc->createAndAddNumbering( seq => $newSeq,
				     seqref => $reference_assmbl_id,
				     refnum => $tile[0],
				     ascending => $ascending );

    # add a file link to the assembly doc in the bsml repository

    $bsmlDoc->createAndAddSeqDataImportN( seq => $newSeq,
				      format => "BSML",
                                      source => $repositoryPath."/".$assmbl_id,
                                      id => "_$assmbl_id");
}

$bsmlDoc->write( $options{'outFile'} );
