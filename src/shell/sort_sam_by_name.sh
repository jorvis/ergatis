#!/bin/sh


INPUT=$1

OUTPUT_DIR=$2

OUTPUT=$3

grep -v "^@" $INPUT > $OUTPUT.reads_only

grep "^@" $INPUT > $OUTPUT.header_only

sort -T $OUTPUT_DIR -k1,1 $OUTPUT.reads_only > $OUTPUT.reads_only.sorted

cat $OUTPUT.header_only $OUTPUT.reads_only.sorted > $OUTPUT.reads.sorted.sam
