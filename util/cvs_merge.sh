#!/bin/bash

# $Id$

# Copyright (c) 2002 The Institute for Genomic Research. All Rights Reserved.

# This script merges a cvs branch back in the trunk

while getopts r:d:b:t opt 2>/dev/null
do case "$opt" in
    r)    RELEASE_TAG=$OPTARG;;
    b)    BRANCH_TAG=$OPTARG;;
    d)    DIRECTORY=$OPTARG;;
    t)    TEST="-n";;
    \?)   echo "Usage: `basename $0` -r(oot tag) -b(branch tag) [-d(irectory for output)] [-t(est)]";
          echo;
          echo "eg. `basename $0` -r papyrus-v1r8 -b papyrus-v1r8-branch";
          exit 1;;
  esac
done

shift `expr $OPTIND - 1`

if [ -z $BRANCH_TAG ]; then
    echo "You must specify a branch tag with -b"
    exit 1;
fi

if [ -z $RELEASE_TAG ]; then
    echo "You must specify a release tag with -r"
    exit 1;
fi

if [ -z $DIRECTORY ]; then
#    echo "Using /tmp as output directory"
    DIRECTORY=/tmp
fi


cd $DIRECTORY
if [ -d merge ]; then
    rm -rf merge
fi
mkdir merge
if [ "$?" -ne 0 ]; then
    echo "Unable to create merge directory"
    exit 1
fi
cd merge
mkdir trunk
mkdir branch
cd trunk
cvs -Q co papyrus_all

cd ../branch
cvs -Q co -r $BRANCH_TAG papyrus_all
if [ "$?" -ne 0 ]; then
    echo "Bad release tag $BRANCH_TAG for papyrus_all"
    exit 1
fi

#Tag the trunk
cd ../trunk
cvs -Q $TEST tag $RELEASE_TAG-premerge-1
cd ../
echo "Added tag $RELEASE_TAG-premerge-1 to trunk"

#Tag the branch
cd branch
cvs -Q $TEST tag $BRANCH_TAG-premerge-1
cd ../
echo "Added $BRANCH_TAG-premerge-1 to branch"

#Merge the branch into the trunk
cd trunk
cvs -q update -kk -j $BRANCH_TAG >& cvsmerge.out
cvs -q status | grep conflict > cvsmerge.conflicts
echo "Merged $BRANCH_TAG into the trunk"

echo "**Inspect $DIRECTORY/merge/trunk/cvsmerge.out and $DIRECTORY/merge/trunk/cvsmerge.conflicts for errors"
echo 
echo "**Commit $DIRECTORY/merge/trunk/"
echo 
echo "**Tag merged code in $DIRECTORY/merge/trunk/: cvs tag $RELEASE_TAG-merge-1"




