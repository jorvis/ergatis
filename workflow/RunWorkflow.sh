#!/bin/sh

umask 0

INSTALL_DIR=/usr/local/devel/ANNOTATION/workflow-sge

while getopts i:l:o:v:h opt
do case "$opt" in
      i) instance=$OPTARG;;
      l) log=$OPTARG;;
      o) out=$OPTARG;;
      v) verbose=$OPTARG;;
      h) echo "Usage: `basename $0` -i instance file -l log file - o output file for stdin and stdout"
         echo;
	 exit;;
      esac
done
	  

source $INSTALL_DIR/exec_env.bash

if [ $verbose ]
then
    debugstr="-v $verbose"
fi

echo "Running workflow instance $instance"

exec RunWorkflow -i $instance -l $log $debugstr --logconf=/usr/local/devel/ANNOTATION/ard/workflow.log4j.properties

echo "To launch viewer execute... WorkflowMonitor $instance"

