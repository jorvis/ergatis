#!/bin/bash

# $Id$

# Copyright (c) 2002 The Institute for Genomic Research. All Rights Reserved.

# This script generates a release tag in cvs
while getopts r:d:t opt 2>/dev/null
do case "$opt" in
    r)    RELEASE_TAG=$OPTARG;;
    d)    DIRECTORY=$OPTARG;;
    t)    TEST="-n";;
    \?)   echo "Usage: `basename $0` -r(release tag) [-d(irectory for output)] [-t(est)]";
          echo;
          echo "eg. `basename $0` -r cram-v1r8b1";
          exit 1;;
  esac
done

if [ -z $RELEASE_TAG ]; then
    echo "You must specify a release tag with -r"
    exit 1;
fi

if [ -z $DIRECTORY ]; then
#    echo "Using /tmp as output directory"
    DIRECTORY=/tmp
fi


cd $DIRECTORY
if [ -d bug_release ]; then
    rm -rf bug_release
fi
mkdir bug_release
if [ "$?" -ne 0 ]; then
    echo "Unable to create bug_release directory"
    exit 1
fi
cd bug_release
cvs -Q co cram_all
cvs -Q $TEST tag $RELEASE_TAG
echo "Added tag $RELEASE_TAG in cram_all"
