#!/bin/bash

# $Id$

# Copyright (c) 2002 The Institute for Genomic Research. All Rights Reserved.

# This script creates a cvs code branch

while getopts d:b:t opt 2>/dev/null
do case "$opt" in
    b)    BRANCH_TAG=$OPTARG;;
    d)    DIRECTORY=$OPTARG;;
    t)    TEST="-n";;
    \?)   echo "Usage: `basename $0` -b(branch tag) [-d(irectory for output)] [-t(est)]";
          echo;
          echo "eg. `basename $0` -b papyrus-v1r8-branch";
          exit 1;;
  esac
done

if [ -z $BRANCH_TAG ]; then
    echo "You must specify a branch tag with -b"
    exit 1;
fi

if [ -z $DIRECTORY ]; then
    #echo "Using /tmp as output directory"
    DIRECTORY=/tmp
fi


cd $DIRECTORY
if [ -d branch ]; then
    rm -rf branch
fi
mkdir branch
if [ "$?" -ne 0 ]; then
    echo "Unable to create branch directory"
    exit 1
fi
cd branch
cvs -Q co papyrus_all
cvs -Q $TEST tag -b $BRANCH_TAG
echo "Created branch $BRANCH_TAG in papyrus_all"


