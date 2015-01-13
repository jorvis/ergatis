#!/bin/sh
# simple_merge_lists.sh- This will merge comma-separated list files together using `cat`
# This is not to be confused with merge_lists.pl, which takes a directory of lists as input
# Author: Shaun Adkins

if [ $# -lt 1 ]
then 
    printf "%b" "Error: No list files provided.\n" >&2
fi
    
while getopts "i:o:" opt
do
	case $opt in
		i) list_string=$OPTARG;;
		o) output=$OPTARG;;
	esac
done    
    
# Take comma-separated string of list files and replace all commas with spaces
temp_string=$(echo $list_string | sed "s/,/ /g")

#Open every list file and cat contents into output file
cat $temp_string > $output
