#!/usr/local/packages/perl-5.8.5/bin/perl

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
use BSML::BsmlBuilder;
}

my %options = ();
GetOptions( \%options, 
	    'tilingPath|t=s', 
	    'referencePath|r=s',
	    'queryPath|q=s',
	    'outFile|o=s',
	    'reference_class=s',
	    'query_class=s',
	    'analysis_id|a=s',
	    'help|h' 
	    ) || pod2usage();

# display documentation
if( $options{'help'} ){
    pod2usage( {-exitval=>0, -verbose => 2, -output => \*STDOUT} );
}

if( !( $options{'tilingPath'} ) ) {
    pod2usage( {-exitval=>1, -verbose => 2, -output => \*STDOUT} );
}

#path to original reference and query mfsa required for SeqDataImport in bsml
if( !( $options{'referencePath'} ) ) {
    pod2usage( {-exitval=>1, -verbose => 2, -output => \*STDOUT} );
}
if( !( $options{'queryPath'} ) ) {
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

## handle some defaults
if (! $options{analysis_id}) {
    $options{analysis_id} = 'unknown_analysis';
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

my $doc = new BSML::BsmlBuilder;
my $ref_id = "";
my $ref_length = -1;
my $reference_class = $options{'reference_class'};
my $query_class = $options{'query_class'};
my %unique_contigs; #track that each contig is only tiled once, otherwise it's ambiguous

## add the analysis element
#$doc->createAndAddAnalysis(
#                            id => $options{analysis_id},
#                            sourcename => $options{'outFile'},
#                          );
#predefine sequences w/ analysis element to track unused sequences
predefineSequences($doc, $options{'referencePath'}, $reference_class);   #reference
predefineSequences($doc, $options{'queryPath'}, $query_class);#query

open( TILINGS, $options{'tilingPath'} ) or die "Unable to open $options{'tilingPath'}";
while (my $line = <TILINGS>) {
    chomp $line;
    #Add in a reference sequence
    if ($line =~ />([\S]*)[\s]([\d]*)/) {
	$ref_id = $1;
	$ref_length = $2;

	if( !( $doc->returnBsmlSequenceByIDR( $ref_id )) ){
	    die "Reference sequence $ref_id missing from predefined sequences";
	}

    }
    #Add in a tile
    else {
	unless ($ref_id) {
	    die "No valid reference sequence";
	}
	my @tile = split("\t",$line);
	my $tile_start = $tile[0] - 1; #interbase 0 coordinates
	my $tile_end = $tile[1]; #not used
	my $tile_length = $tile[3];
	my $percent_coverage = $tile[4];
	my $percent_identity = $tile[5];
	my $tile_id = $tile[7];

	my $tile_ascending = 1;
	my $tile_is_complement = 0;
	my $refnum = $tile_start;

	if( $tile[6] eq '-' ) { #was on the complement
	    $tile_is_complement = 1;
	    $tile_ascending = 0;
	    $refnum = $tile_end;
	}

	#ensure uniqueness
	if (exists $unique_contigs{$tile_id}) {
	    die "Ambigous multi-tiling of contig $tile_id";
	}
	else { $unique_contigs{$tile_id} = 1; }

	#check that tile.length < scaffold.length
	if ($ref_length < $tile_length) {
	    die "Length of scaffold $ref_id ($ref_length), less than tile ($tile_id) ($tile_length)";
	}

	#negative coords are okay, they occur with -c and tiling across origin
	#if ($tile_start < 0 || $tile_length < 0) {
	#    die "Invalid tile_start ($tile_start) or tile_length ($tile_length) for $tile_id in $options{'tilingPath'}";
	#}

	if( !( $doc->returnBsmlSequenceByIDR( $tile_id )) ){
	    die "Query sequence $tile_id missing from predefined sequences";
	}
	else {
	    my $tile_sequence = $doc->returnBsmlSequenceByIDR( $tile_id );

	    #Numbering element only applicable for exact matches
	    if ($percent_identity == 100 && $percent_coverage == 100) {		
		$doc->createAndAddNumbering( seq => $tile_sequence,
					     seqref => $ref_id,
					     refnum => $refnum,
					     ascending => $tile_ascending );
	    }

	    my $aln = $doc->createAndAddSequencePairAlignment(
                                         refseq => $ref_id,
					 compseq => $tile_id,
					 method=> 'tiling',
                                         );
	    $aln->addBsmlLink('analysis', '#' . $options{analysis_id}, 'computed_by');

	    my $run = $doc->createAndAddSequencePairRun(
                                           alignment_pair => $aln,
                                           refpos => $tile_start,
                                           runlength => $tile_length,
                                           refcomplement => 0, #ref sequence is never complement
					   #tile is always entire contig
                                           comppos => 0,
                                           #comprunlength => $tile_sequence->returnattr( 'length' ),
                                           comprunlength => $tile_length,
                                           compcomplement => $tile_is_complement );

	    #add percent_identity and percent_coverage Attributes
	    $run->addBsmlAttr( 'percent_identity', $percent_identity);
	    $run->addBsmlAttr( 'percent_coverage', $percent_coverage);
	}	
    }
}
close(TILINGS);

$doc->write( $options{'outFile'} );


sub predefineSequences {
    my $doc = shift;
    my $seqfile = shift;
    my $class = shift;

    my %seqs = loadMultiSequence( $seqfile );

    for my $seqid ( sort {$a<=>$b} keys %seqs ) {
        ## capture the first element of the header up to the first whitespace
        my $id;
        if ($seqs{$seqid}{h} =~ /^(\S+)/) {
            $id = $1;
        } else {
            die("unrecognized header format: $seqs{$seqid}{h}");
        }

	my $seqstub = $doc->createAndAddSequence(
						 $id, #id
						 undef, #title
						 length($seqs{$seqid}{s}), #length
						 'dna', #molecule
						 $class #class
						 );
	#Seq-data-import/@identifier must equal the fasta header up to the first space	    
	my $dataimport = $doc->createAndAddSeqDataImport(
 				    $seqstub,              # Sequence element object reference
 				    'fasta',               # //Seq-data-import/@format
 				    $seqfile,              # //Seq-data-import/@source
 				    undef,                 # //Seq-data-import/@id
 				    $id                    # //Seq-data-import/@identifier
 				    );
	$seqstub->addBsmlLink('analysis', '#' . $options{analysis_id}, 'input_of');
    }
}

#taken from fasta2bsml.pl
sub loadMultiSequence {
    #  USAGE:   loadMultiSequence($filepath)
    #  RETURNS: hash
    #
    #  takes a file or path as an argument.  that file should be a multiple-
    #  sequence FASTA file.  It returns a hash with a structure like:
    #      $db{id}{'h'} = header
    #             {'s'} = sequence without whitespace
    #
    #  where id is an incrementing integer that represents that sequence's
    #  order in the file.
    #
    #########################################################################
    my ($file) = @_;
    
    my $seqid = 0;
    my $seq = '';
    my $header;
    my %db;
    
    ## load the sequence file
#    open (my $sfh, "<$file") || $logger->logdie("can't open $file because $!");
    open (my $sfh, "<$file") || die("can't open $file because $!");

    for (<$sfh>) {
        ## if we find a header line ...
        if (/^\>(.*)/) {

            $header = $1;

            ## don't do anything if this is the first sequence
            if ($seqid == 0) {
                $seqid++;
                $db{$seqid}{'h'} = $header;
                next;
            } 

            ## remove whitespace
            $seq =~ s/\s//g;
 
            ## record the previous sequence before starting the new one
            $db{$seqid}{'s'} = $seq;

            ## increment the id counter
            $seqid++;

            ## record the new header
            $db{$seqid}{'h'} = $header;

            ## reset the sequence
            $seq = '';

        ## else we've found a sequence line
        } else {
            ## skip it if it is just whitespace
            next if (/^\s*$/);

            ## record this portion of the sequence
            $seq .= $_;
        }
    }
    
    ## don't forget the last sequence
    $seq =~ s/\s//g;
    $db{$seqid}{'s'} = $seq;

    ## close the sequence file
    close $sfh;
    
    return %db;
}
