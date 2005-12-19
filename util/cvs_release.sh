#!/bin/bash

# $Id$

# Copyright (c) 2002 The Institute for Genomic Research. All Rights Reserved.

# This script generates a release tag in cvs
while getopts r:d:pt opt 2>/dev/null
do case "$opt" in
    r)    RELEASE_TAG=$OPTARG;;
    d)    DIRECTORY=$OPTARG;;
    p)    PRINT="yes";;
    t)    TEST="-n";;
    \?)   echo "Usage: `basename $0` -r(release tag) [-d(irectory for output)] [-t(est)] [-p(rint latest release]";
          echo;
          echo "eg. `basename $0` -r ergatis-v1r8b1";
          exit 1;;
  esac
done

if [ $PRINT ]; then
    cvs log cvs_release.sh | perl -e 'while(my $line=<STDIN>){if($line=~/symbolic/){$retrieve=1;}if($retrieve==1){($tag) = ($line =~ /(ergatis-v\d+r\d+b\d+)\:/);if($tag){print $tag,"\n";exit}}}'
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
if [ -d bug_release ]; then
    rm -rf bug_release
fi
mkdir bug_release
if [ "$?" -ne 0 ]; then
    echo "Unable to create bug_release directory"
    exit 1
fi
cd bug_release
cvs -Q co ergatis_all
cvs -Q $TEST tag $RELEASE_TAG
echo "Added tag $RELEASE_TAG in ergatis_all"
