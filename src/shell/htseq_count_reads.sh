#!/bin/sh

M=$1

T=$2

A=$3

I=$4

S=$5

GTF=$6

SAM=$7

COUNTFILE=$8

export PATH=/usr/local/packages/Python-2.6.4/bin:$PATH
export PYTHONPATH=/usr/local/packages/pythonlib/2.6

#use python2.6

python -m HTSeq.scripts.count -m $M -t $T -i $I -s $S $SAM $GTF > $COUNTFILE