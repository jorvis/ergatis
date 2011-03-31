#!/bin/sh

sbin=$1
shift
odir=$1
shift
alignment_file=$1
shift
paired_end_file=$1
shift
sizes_file=$1
shift
chrom_seq_dir=$1
shift
min_splice_support=$1
shift

echo "$sbin $odir $alignment_file $paired_end_file $sizes_file $chrom_seq_dir $min_splice_support"

other_args=""
while(( "$#" )); do
other_args=$other_args" "$1
shift
done

for i in `awk '{print $1}' $sizes_file` 
do 
mkdir $odir/$i
cd $odir/$i
java -Xmx4000m -jar $sbin -alignment $alignment_file -pairedEnd $paired_end_file -sizeFile $sizes_file -chr \'$i\' -chrSequence $chrom_seq_dir/$i.fa -out $i.segments -minSpliceSupport $min_splice_support $other_args
done

