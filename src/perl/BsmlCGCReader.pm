package BsmlCGCReader;
@ISA = qw( BSML::BsmlReader );

use strict;
use BSML::BsmlReader;

# This workflow "helper" function may better serve in the client layer as it is really just
# a summary of the coordinates produced by geneIdtoGenomicCoords.

# For each gene on the input assembly, returns the coordinates of the gene, 
# namely the genomic positions that span the gene's largest transcript. Results
# are returned as a reference to anonymous hashes.  

sub fetch_gene_positions
  {
    my $self = shift;
    my ($assemblyId) = @_;

    my $rlist = [];

    foreach my $gene ( @{$self->assemblyIdtoGeneList($assemblyId)} )
      {
	my $rhash = {};
        my $geneCoords = $self->geneIdtoGenomicCoords( $gene );

	$rhash->{$gene}->{'startpos'} = $geneCoords->[0]->{'GeneSpan'}->{'startpos'};
	$rhash->{$gene}->{'endpos'} = $geneCoords->[0]->{'GeneSpan'}->{'endpos'};
	$rhash->{$gene}->{'complement'} = $geneCoords->[0]->{'GeneSpan'}->{'complement'};

	push( @{$rlist}, $rhash );
      }

    return $rlist;
  }

#this function might serve better as a bitscore uitility function than in the general API

sub fetchAlignmentScoresBetweenAssemblies {

    my $reader = shift;
    my $query_asmbl_id = shift;
    my $match_asmbl_id = shift;
    #my $order = shift || '1';

    my $bit_score_hash={};
    foreach my $seq (@{$reader->assemblyIdtoSeqList($query_asmbl_id )}) {   #grab all seq obj given query_asmbl_id
	my $seq_id = $seq->returnattr( 'id' );                              #return seq_id of an seq obj
	foreach my $aln (@{$reader->fetch_all_alignmentPairs( $seq_id )}) { #return all alignment with query as $seq_id
	    if(ref($aln)) {
		my $match_ref = $reader->readSeqPairAlignment($aln);            #return all pair_runs for an alignment_pair
		my $m_asmbl_id = $reader->seqIdtoAssemblyId($match_ref->{'compseq'});
		if($m_asmbl_id eq $match_asmbl_id) {                            #check to see if match gene belongs to match asmbl_id
		    my $best_bit_score=0;
		    foreach my $pair_run(@{ $match_ref->{'seqPairRuns'} }) {
			$best_bit_score = $pair_run->{'runscore'} if($pair_run->{'runscore'} > $best_bit_score);  #store best bit_score
		    }
		    $bit_score_hash->{$seq_id}->{ $match_ref->{'compseq'} }->{'bit_score'} = $best_bit_score;
		  }
		#add the best bit score for a query gene against itself 
		#to provide a baseline for bit score comparison
		#-----------------------------------------------------
		my $aln = $reader->fetch_all_alignmentPairs($seq_id, $seq_id);
		$match_ref = $reader->readSeqPairAlignment($aln->[0]);            #return all pair_runs for an alignment_pair
		my $best_bit_score=0;
		foreach my $pair_run(@{ $match_ref->{'seqPairRuns'} }) {
		    $best_bit_score = $pair_run->{'runscore'} if($pair_run->{'runscore'} > $best_bit_score);  #store best bit_score
		}
		$bit_score_hash->{$seq_id}->{$seq_id}->{'bit_score'} = $best_bit_score;
		#------------------------------------------------------
	    }
	}

    }

	return $bit_score_hash;
}

sub fetch_genome_pairwise_matches
  {
    my $reader = shift;
    my ($query_asmbl_id, $match_asmbl_id) = @_;

    my $rlist=[];

    foreach my $seq (@{$reader->assemblyIdtoSeqList($query_asmbl_id )}) {   #grab all seq obj given query_asmbl_id
	my $seq_id = $seq->returnattr( 'id' );                              #return seq_id of an seq obj
	foreach my $aln (@{$reader->fetch_all_alignmentPairs( $seq_id )}) { #return all alignment with query as $seq_id
	    if(ref($aln)) {
		my $match_ref = $reader->readSeqPairAlignment($aln);            #return all pair_runs for an alignment_pair
		my $m_asmbl_id = $reader->seqIdtoAssemblyId($match_ref->{'compseq'});
		my $q_asmbl_id = $reader->seqIdtoAssemblyId($match_ref->{'refseq'});
		if($m_asmbl_id eq $match_asmbl_id || ($match_asmbl_id eq 'all' && ($m_asmbl_id ne $q_asmbl_id) ))
		  {
		    my $rhash = {};
		    
		    $rhash->{'query_gene_name'} = $match_ref->{'refseq'};
		    $rhash->{'match_gene_name'} = $match_ref->{'compseq'};

		    my $best_percent_identity = 0.0;
	      
		    foreach my $pair_run(@{ $match_ref->{'seqPairRuns'} }) {


			$best_percent_identity = $pair_run->{'percent_identity'} if($pair_run->{'percent_identity'} > $best_percent_identity);
		      }

		    $rhash->{'percent_identity'} = $best_percent_identity;

		    my $best_percent_similarity = 0.0;
		    
		    foreach my $pair_run(@{ $match_ref->{'seqPairRuns'} }) {
			$best_percent_similarity = $pair_run->{'percent_similarity'} if($pair_run->{'percent_similarity'} > $best_percent_similarity);
		      }

		    $rhash->{'percent_similarity'} = $best_percent_similarity;

		    my $best_pval = 1e30;
		    foreach my $pair_run(@{ $match_ref->{'seqPairRuns'} }) {
		      $best_pval = $pair_run->{'p_value'} if($pair_run->{'p_value'} < $best_pval);
		    }

		    $rhash->{'pval'} = $best_pval;

		    push( @{$rlist}, $rhash );
		  }
		  
	      }
	  }
      }
    return $rlist;
  }

1;
