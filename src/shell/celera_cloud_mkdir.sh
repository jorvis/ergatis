#!/bin/sh

dir=$1

if [ -x $dir ] 
then
	echo "Directory $dir already exists, will not create it"
	exit 0
else
	mkdir $dir
fi

