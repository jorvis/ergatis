#!/bin/sh

INSTALL_DIR=/usr/local/devel/ANNOTATION/wf-test

while getopts t:c:i:l:o:v:h opt
do case "$opt" in
      t) template=$OPTARG;;
      c) ini=$OPTARG;;
      i) instance=$OPTARG;;
      l) log=$OPTARG;;
      o) out=$OPTARG;;
      v) verbose=$OPTARG;;
      h) echo "Usage: `basename $0` -t template file -c ini file -i instance file -l log file - o output file for stdin and stdout"
         echo;
	 exit;;
      esac
done
	  

source $INSTALL_DIR/exec_env.bash

if [ $verbose ]
then
    debugstr="-v $verbose"
fi

echo "Creating workflow instance $instance"
CreateWorkflow -t $template -c $ini -i $instance -l $log $debugstr --autobuild=false 1> $out 2> $out.stderr


