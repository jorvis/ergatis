#!/bin/sh

while getopts d: opt
do case "$opt" in
  d) directory=$OPTARG;;
  esac;
done
ps | grep java | perl -ne '($id) = ($_ =~ /\s*(\d+)/);`kill -2 $id`;'

while ps | grep java 
do
  echo .
done

if [ -d $directory ]
then
    find $directory -name "*.xml" -exec perl -pi -e 's/>running/>interrupted/' \{\} \;
    find $directory -name "*.xml" -exec perl -pi -e 's/>pending/>interrupted/' \{\} \;
fi

