#! /bin/bash


# These are the parameters for find_intergenic_background_cutoff.py.  I've
# explicitly specified the default values of the three optional parameters
# just to show what the option names are
quantile=0.7
min_interbutt=50
max_interbutt=1000

# These are the remaining parameters of find_transcripts.pl
min_transcript=100
base_feature_type=CDS
max_intron=0
min_intergenic=1
stranded=1

###############################################################################
# command line parameter

bin_d=$1
shift
output_d=$1
shift
wig_d=$1
shift
reference_fasta_f=$1
shift
gff_f=$1
shift
sample=$1
shift
QUAN=${1:-$quantile}
shift
MIN_IB=${1:-$min_interbutt}
shift
MAX_IB=${1:-$max_interbutt}
shift
MIN_TS=${1:-$min_transcript}
shift
FEATURE=${1:-$base_feature_type}
shift
MAX_IN=${1:-$max_intron}
shift
MIN_IG=${1:-$min_intergenic}
shift
STRAND=${1:-$stranded}
shift

###############################################################################

contig_sizes_f=$output_d/contig_sizes.tsv
< $reference_fasta_f /home/edrabek/bin/fastalen.awk \
  > $contig_sizes_f

$bin_d/find_intergenic_background_cutoff.py \
    --quantile=$QUAN \
    --min_interbutt=$MIN_IB \
    --max_interbutt=$MAX_IB \
    $contig_sizes_f \
    $gff_f \
    $wig_d/$sample.{forward,reverse}_coverage.wig \
    > $output_d/$sample.cutoff

cutoff=`cat $output_d/$sample.cutoff`

perl $bin_d/find_transcripts.pl \
    --file_type=WIG \
    --file=$wig_d/$sample.forward_coverage.wig,$wig_d/$sample.reverse_coverage.wig \
    --sizes_file=$contig_sizes_f \
    --stranded=$STRAND \
    --min_cov=$cutoff \
    --max_intron=$MAX_IN \
    --min_transcript=$MIN_TS \
    --min_intergenic=$MIN_IG \
    --base_gff3=$gff_f \
    --base_feature_type=$FEATURE \
    --output_stub=$output_d/$sample.transcripts
