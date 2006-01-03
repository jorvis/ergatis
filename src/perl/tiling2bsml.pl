#!/usr/local/packages/perl-5.8.5/bin/perl

eval 'exec /usr/local/packages/perl-5.8.5/bin/perl  -S $0 ${1+"$@"}'
    if 0; # not running under some shell

eval 'exec /usr/local/packages/perl-5.8.5/bin/perl  -S $0 ${1+"$@"}'
    if 0; # not running under some shell

# this is temporary for development
# use lib "/export/CVS/bsml/src/";

=head1  NAME 

tiling2bsml.pl  - convert show_tiling output into BSML sequence mapping documents

=head1 SYNOPSIS

USAGE:  tiling2bsml.pl -t tilings_file -b path_to_bsml_repository -o output_file_bsml

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

tiling2bsml.pl is designed to convert tiling path data from show_tiling 
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
	    'reference_class=s',
	    'query_class=s',
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

#needed inputs available?
if( !( $options{'reference_class'} ) ) {
    pod2usage( {-exitval=>1, -verbose => 2, -output => \*STDOUT} );
}
if( !( $options{'query_class'} ) ) {
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
my $ref_length = -1;
my $reference_class = $options{'reference_class'};
my $query_class = $options{'query_class'};
my %unique_contigs; #track that each contig is only tiled once, otherwise it's ambiguous
while (my $line = <TILINGS>) {
    chomp $line;
    #Add in a reference sequence
    if ($line =~ />([\S]*)[\s]([\d]*)/) {
	$ref_id = $1;
	$ref_length = $2;

	my $ref_sequence = $doc->createAndAddSequence(
				    $ref_id, #id
				    undef, #title
				    $ref_length, #length
				    'dna', #molecule
				    $reference_class #class
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
	if (exists $unique_contigs{$tile_id}) {
	    die "Ambigous multi-tiling of contig $tile_id";
	}
	else { $unique_contigs{$tile_id} = 1; }

	#check that tile.length < scaffold.length
	if ($ref_length < $tile_length) {
	    die "Length of scaffold $ref_id ($ref_length), less than tile ($tile_id) ($tile_length)";
	}

	#check that the coordinates are positive
	if ($tile_start < 0 || $tile_length < 0) {
	    die "Invalid tile_start ($tile_start) or tile_length ($tile_length) for $tile_id in $options{'tilingPath'}";
	}
	

	my $tile_sequence = $doc->createAndAddSequence(
				    $tile_id, #id
				    undef, #title
				    $tile_length, #length
				    'dna', #molecule
				    $query_class #class
				    );
	#Seq-data-import/@identifier must equal the fasta header up to the first space	    
	my $tile_sequence_data = $doc->createAndAddSeqDataImport(
 				    $tile_sequence,              # Sequence element object reference
 				    'fasta',                    # //Seq-data-import/@format
 				    $options{'tilingPath'},               # //Seq-data-import/@source
 				    undef,                      # //Seq-data-import/@id
 				    $tile_id                     # //Seq-data-import/@identifier
 				    );
	$doc->createAndAddNumbering( seq => $tile_sequence,
 				     seqref => $ref_id,
 				     refnum => $tile_start,
 				     ascending => $tile_ascending );	
    }
}
close(TILINGS);

$doc->write( $options{'outFile'} );
