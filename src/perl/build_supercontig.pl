#!/usr/local/packages/perl-5.8.5/bin/perl
#run as:
# perl -I /usr/local/devel/ANNOTATION/ard/chado-v1r5b1/lib/site_perl/5.8.5/ build_supercontig.pl --bsml_input=/usr/local/scratch/annotation/AGUSSMAN/output_repository/tiling/4563_50-65/tiling.bsml --bsml_output=./test.bsml --scaffold_class=ultracontig --fasta_output=./test.fsa
#or
#perl -I /usr/local/devel/ANNOTATION/ard/chado-v1r5b1/lib/site_perl/5.8.5/ bin/build_supercontig.pl --bsml_input=/tmp/agussman/tiling.bsml --bsml_output=/tmp/agussman/build.bsml --scaffold_class=ultracontig --fasta_output=/tmp/agussman/build.fsa

eval 'exec /usr/local/packages/perl-5.8.5/bin/perl  -S $0 ${1+"$@"}'
    if 0; # not running under some shell

#
# perldoc goes here
#
=head1  NAME 

build_supercontig.pl - convert ordered Seq-pair-alignments into a single contig (most likely to be run on the output of tiling2bsml)
=cut

use strict;
use warnings;

use Getopt::Long qw(:config no_ignore_case no_auto_abbrev);
use Pod::Usage;
umask(0000);


BEGIN {
use BSML::BsmlBuilder;
use Workflow::Logger;
use BSML::BsmlReader;
use BSML::BsmlParserSerialSearch;
use BSML::BsmlParserTwig;
    
}

my %options = ();
GetOptions( \%options, 
	    'bsml_input|i=s',
	    'bsml_output|o=s',
	    'fasta_output=s',
	    'analysis_id|a=s',
	    'scaffold_class=s',
	    'help|h' 
	    ) || pod2usage();

# display documentation
if( $options{'help'} ){
    pod2usage( {-exitval=>0, -verbose => 2, -output => \*STDOUT} );
}
if( !( $options{'bsml_input'} ) ) {
    pod2usage( {-exitval=>1, -verbose => 2, -output => \*STDOUT} );
}
if( !( $options{'bsml_output'} ) ) {
    pod2usage( {-exitval=>1, -verbose => 2, -output => \*STDOUT} );
}
if( !( $options{'fasta_output'} ) ) {
    pod2usage( {-exitval=>1, -verbose => 2, -output => \*STDOUT} );
}

#handle some defaults
if (! $options{analysis_id}) {
    $options{analysis_id} = 'unknown_analysis';
}
if (! $options{scaffold_class}) {
    $options{scaffold_class} = 'supercontig';
}

#create output file
my $fasta_out = $options{'fasta_output'};
open (my $FFILE, ">$fasta_out") || die "Unable to open for writing $fasta_out";
close($FFILE);

#the spacer sequence following each tile in supercontig
my $spacer_seq = "NNNNNCACACACTTAATTAATTAAGTGTGTGNNNNN";
my $spacer_len = length($spacer_seq);

my $bsmlFile = $options{'bsml_input'};
my %tiles;
my %seqtrack; #track sequence elements for removal of sequence data later
my $doc = new BSML::BsmlBuilder;
# add the analysis element stub
$doc->createAndAddAnalysis(
			   id => $options{analysis_id},
			   sourcename => $options{bsml_output},
			   );

#obtain alignment info (for building supercontigs) from input bsml
#obtain sequence info from input bsml
my $alnParser = new BSML::BsmlParserSerialSearch( AlignmentCallBack => \&alignmentHandler );
my $seqParser = new BSML::BsmlParserSerialSearch( SequenceCallBack => \&sequenceHandler );
$seqParser->parse( $bsmlFile );
$alnParser->parse( $bsmlFile );

#build a supercontig for each reference scaffold
my $reader = new BSML::BsmlReader();
foreach my $refseq (keys %tiles) {
    my $supercontig_seq = "";
    my $supercontig_name = $options{scaffold_class}."_".$refseq;
    my $supercontig_len = 0; #current length of supercontig
    #add stub supercontig sequence because of catch22 in adding alignments
    my $seqstub = $doc->createAndAddSequence(
 					     $supercontig_name, #id
 					     undef, #title
 					     0, #length
 					     'dna', #molecule
 					     $options{scaffold_class} #class
 					     );
    $seqstub->addBsmlLink('analysis', '#' . $options{analysis_id}, 'computed_by');

    #build supercontig in order of:
    # 1. start position on reference scaffold
    # 2. length of tile
    # in actuality, only 1 is relevant as show-tilings doesn't output subsumed sequences
    foreach my $refpos (sort { $a <=> $b } keys %{$tiles{$refseq}}) {
	foreach my $runlength (sort {$a<=>$b} keys %{$tiles{$refseq}{$refpos}}) {
	    #foreach my $compseq (@{$tiles{$refseq}{$refpos}{$runlength}}) {
	    foreach my $compseq (keys %{$tiles{$refseq}{$refpos}{$runlength}}) {
		my $is_comp = $tiles{$refseq}{$refpos}{$runlength}{$compseq};
		my $refnum = $supercontig_len;
		my $is_ascending = 1;
		print "$refseq\t$refpos\t$runlength\t$compseq\t$is_comp\t($supercontig_name)\n";

		#add to supercontig sequence
		my $seq = $reader->subSequence(${$seqtrack{$compseq}},-1,0,0);
		$seq =~ s/\s//g;
		my $seq_len = length($seq);
		if ($seq_len < 1) {
		    die "Grevious error regarding the recommended length of a sequence";
		}

		#check if on complement (Numbering element is at end, descending; Sequence is revcomp)
		if ($is_comp) {
		    $is_ascending = 0;
		    $refnum += $seq_len;
		    $seq = &revcomp($seq);
		}

		$seq .= $spacer_seq; #add spacer

		#add Numbering
		my $seq_elm = $doc->returnBsmlSequenceByIDR( $compseq ) or die "Unable to retrieve sequence";
		$doc->createAndAddNumbering( seq => $seq_elm,
					     seqref => $supercontig_name,
					     refnum => $refnum,
					     ascending => $is_ascending );

		#add SeqPairAlignment
		my $aln = $doc->createAndAddSequencePairAlignment(
								  refseq => $supercontig_name,
								  compseq => $compseq,
								  method=> 'tiling',
								  );
		$aln->addBsmlLink('analysis', '#' . $options{analysis_id}, 'computed_by');
 		my $run = $doc->createAndAddSequencePairRun(
                                            alignment_pair => $aln,
                                            refpos => $supercontig_len,
                                            runlength => $seq_len,
                                            refcomplement => 0, #ref sequence is never complement
                                            comppos => 0,       #always entire contig
                                            comprunlength => $seq_len,
					    compcomplement => $is_comp );
		#by definition these are 100%
		$run->addBsmlAttr( 'percent_identity', '100.00');
		$run->addBsmlAttr( 'percent_coverage', '100.00');

		$supercontig_len += $seq_len + $spacer_len;
		$supercontig_seq .= $seq;
	    }
	}
    }

    #output supercontig to file and update the supercontig Sequence element
    #append output file
    open (my $FFILE, ">>$fasta_out") || die "Unable to append to $fasta_out";
    print {$FFILE} fasta_out($supercontig_name, $supercontig_seq);
    close($FFILE);

    $seqstub->setattr( 'length', $supercontig_len);
    my $dataimport = $doc->createAndAddSeqDataImport(
 				    $seqstub,              # Sequence element object reference
 				    'fasta',               # //Seq-data-import/@format
 				    $fasta_out,            # //Seq-data-import/@source
 				    undef,                 # //Seq-data-import/@id
 				    $supercontig_name      # //Seq-data-import/@identifier
 				    );
}

#output results
$doc->write( $options{'bsml_output'} );
chmod(0666, $options{'bsml_output'});
chmod(0666, $fasta_out);

#
# subroutines follow
#

# Create stubs for all the sequences in the input bsml
sub sequenceHandler {
    my $oldseq = shift;
    print "Found sequence!";
    #pull attributes from original sequence
    my $id = $oldseq->returnattr( 'id' );
    my $length = $oldseq->returnattr( 'length' );
    my $molecule = $oldseq->returnattr( 'molecule' );
    my $class = $oldseq->returnattr( 'class' );

    #create new sequence stub
    my $seqstub = $doc->createAndAddSequence(
 					     $id, #id
 					     undef, #title
 					     $length, #length
 					     $molecule, #molecule
 					     $class #class
 					     );

    #add a sequence data element
    my $seqdataimport = $oldseq->returnBsmlSeqDataImport(); #always returns true, so test id
    if ($seqdataimport->{id}) {
	print " SeqDataImport";
	my $dataimport = $doc->createAndAddSeqDataImport(
 				    $seqstub,              # Sequence element object reference
 				    $seqdataimport->{format},               # //Seq-data-import/@format
 				    $seqdataimport->{source},              # //Seq-data-import/@source
 				    undef,                 # //Seq-data-import/@id
 				    $seqdataimport->{identifier}                    # //Seq-data-import/@identifier
 				    );
    }  
    elsif ( my $seqdata = $oldseq->returnSeqData() ){
	print " SeqData";
	$seqstub->addBsmlSeqData( $seqdata );
    }
#    else {
#	die "No Seq-data or Seq-data-import found";
#    }
    print "\n";

    $seqtrack{$id} = \$oldseq;
    $seqstub->addBsmlLink('analysis', '#' . $options{analysis_id}, 'input_of');
    
}


#build the tile tracking data structure
#order by {refseq}{startpos}{length} ( TileA, TileB,...)
sub alignmentHandler {
    my $aln = shift;

    print "Found alignment!";
    my $refseq = $aln->returnattr( 'refseq');
    my $compseq = $aln->returnattr( 'compseq');
	
    foreach my $SeqPairRun ( @{$aln->returnBsmlSeqPairRunListR()} ) {
	my $runlength = $SeqPairRun->returnattr( 'runlength' );
	my $refpos = $SeqPairRun->returnattr( 'refpos' );
	#my $comprunlength = $SeqPairRun->returnattr( 'comprunlength');
	#my $refcomplement = $SeqPairRun->returnattr( 'refcomplement' );
	#my $comppos = $SeqPairRun->returnattr( 'comppos' );
	my $compcomplemnt = $SeqPairRun->returnattr( 'compcomplement' );
	#actually, there can only ever be one tile beginning at a refpos
	#push(@{$tiles{$refseq}->{$refpos}->{$runlength}}, $compseq);
	$tiles{$refseq}->{$refpos}->{$runlength}->{$compseq} = $compcomplemnt;
    }

    print " $refseq $compseq ";

    print "\n";

}


sub fasta_out {
    #This subroutine takes a sequence name and its sequence and
    #outputs a correctly formatted single fasta entry (including newlines).
    my ($seq_name, $seq) = @_;
    my $fasta=">"."$seq_name"."\n";
    $seq =~ s/\s+//g;
    for(my $i=0; $i < length($seq); $i+=60){
        my $seq_fragment = substr($seq, $i, 60);
        $fasta .= "$seq_fragment"."\n";
    }
    return $fasta; 
}

sub revcomp{
    my($seq) = @_;
    $seq =~ tr/ATGCatgc/TACGtacg/;
    $seq = reverse($seq);
    return $seq;
}
