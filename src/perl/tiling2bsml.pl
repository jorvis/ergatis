#!/usr/local/packages/perl-5.8.5/bin/perl

eval 'exec /usr/local/packages/perl-5.8.5/bin/perl  -S $0 ${1+"$@"}'
    if 0; # not running under some shell

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
    require '/usr/local/devel/ANNOTATION/cas/lib/site_perl/5.8.5/BSML/BsmlBuilder.pm';
    import BSML::BsmlBuilder;
}

my %options = ();
GetOptions( \%options, 
	    'tilingPath|t=s', 
	    'referencePath|r=s', 
	    'outFile|o=s', 
	    'help|h' 
	    ) || pod2usage();

# display documentation
if( $options{'help'} ){
    pod2usage( {-exitval=>0, -verbose => 2, -output => \*STDOUT} );
}

if( !( $options{'tilingPath'} ) ) {
    pod2usage( {-exitval=>1, -verbose => 2, -output => \*STDOUT} );
}

#path to original reference mfsa required for SeqDataImport in bsml
if( !( $options{'referencePath'} ) ) {
    pod2usage( {-exitval=>1, -verbose => 2, -output => \*STDOUT} );
}

if( !( $options{'outFile'} ) ) {
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

#my $repositoryPath = $options{'bsml_repository_path'};

open( TILINGS, $options{'tilingPath'} ) or die "Unable to open $options{'tilingPath'}";
my $doc = new BSML::BsmlBuilder;
my $ref_id = "";
my %unique_contigs; #track that each contig is only tiled once, otherwise it's ambiguous
while (my $line = <TILINGS>) {
    chomp $line;
    #Add in a reference sequence
    if ($line =~ />([\S]*)[\s]([\d]*)/) {
	$ref_id = $1;
	my $ref_length = $2;

	my $ref_sequence = $doc->createAndAddSequence(
				    $ref_id, #id
				    undef, #title
				    $ref_length, #length
				    'mol-not-set', #molecule
				    'supercontig' #class
				    );
	#Seq-data-import/@identifier must equal the fasta header up to the first space	    
	my $ref_sequence_data = $doc->createAndAddSeqDataImport(
 				    $ref_sequence,              # Sequence element object reference
 				    'fasta',                    # //Seq-data-import/@format
 				    $options{'referencePath'},               # //Seq-data-import/@source
 				    undef,                      # //Seq-data-import/@id
 				    $ref_id                     # //Seq-data-import/@identifier
 				    );
    }
    #Add in a tile
    else {
	unless ($ref_id) {
	    die "No valid reference sequence";
	}
	my @tile = split("\t",$line);
	my $tile_id = $tile[7];
	my $tile_length = $tile[3];
	my $tile_start = $tile[0] - 1; #interbase 0 coordinates
	my $tile_ascending = 0;

	if( $tile[6] eq '+' ) {
	    $tile_ascending = 1;
	}
        #  else stay 0

	#ensure uniqueness
#	(exists $unique_contigs{$tile_id}) ? die "Ambigous multi-tiling of contig $tile_id" : $unique_contigs{$tile_id} = 1;
	if (exists $unique_contigs{$tile_id}) {
	    die "Ambigous multi-tiling of contig $tile_id";
	}
	else { $unique_contigs{$tile_id} = 1; }

	my $tile_sequence = $doc->createAndAddSequence(
				    $tile_id, #id
				    undef, #title
				    $tile_length, #length
				    'mol-not-set', #molecule
				    'contig' #class
				    );
	$doc->createAndAddNumbering( seq => $tile_sequence,
 				     seqref => $ref_id,
 				     refnum => $tile_start,
 				     ascending => $tile_ascending );	
    }
}
close(TILINGS);

$doc->write( $options{'outFile'} );

# add the reference sequence to the Tiling Document

# my $reference_assmbl_id = '';

# my $line = <TILINGS>;

# if( !($line) )
# {
#     die "Tiling path contains no reference sequence\n";
# }
# else
# {
#     # Grab the refseq id from the first line
#     if ($line =~ />([\S]*)[\s]([\d]*)/)
#     {
# 	$reference_assmbl_id = $1;

# 	# add the reference sequence to the bsml tiling
# 	my $newSeq = $bsmlDoc->createAndAddExtendedSequenceN( id => $1,
# 							      title => $1,
# 							      molecule => 'dna',
# 							      length => $2 );

# 	$newSeq->addattr( 'class', 'assembly' );

# 	# add a file link to the assembly bsml doc in the repository

# 	$bsmlDoc->createAndAddSeqDataImportN( seq => $newSeq,
# 				      format => "BSML",
#                                       source => $repositoryPath."/".$reference_assmbl_id.".bsml",
#                                       id => "_$reference_assmbl_id");
	
#     }
#     else
#     {
# 	die "Unable to determine reference sequence\n";
#     }
# }

# my $rank = 0;

# # add each tiled sequence to the Tiling Document

# while( my $line = <TILINGS> )
# {
#     my @tile = split( "\t", $line );
#     chomp( $tile[7] );

#     my $assmbl_id = $tile[7];
#     my $orientation = $tile[6];
#     my $ascending = 0;

#     if( $orientation eq '+' )
#     {
# 	$ascending = 1;
#     }
#     else
#     {
# 	$ascending = 0;
#     }

#     my $newSeq = $bsmlDoc->createAndAddExtendedSequenceN( id => $assmbl_id,
# 							  title => $assmbl_id,
# 							  molecule => 'dna',
# 							  length => $tile[3] );

#     $newSeq->addattr( 'class', 'assembly' );

#     $bsmlDoc->createAndAddNumbering( seq => $newSeq,
# 				     seqref => $reference_assmbl_id,
# 				     refnum => $tile[0],
# 				     ascending => $ascending );

#     # add a file link to the assembly doc in the bsml repository

#     $bsmlDoc->createAndAddSeqDataImportN( seq => $newSeq,
# 				      format => "BSML",
#                                       source => $repositoryPath."/".$assmbl_id,
#                                       id => "_$assmbl_id");
# }

# $bsmlDoc->write( $options{'outFile'} );
