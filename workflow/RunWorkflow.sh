#!/bin/sh

umask 0

INSTALL_DIR=/usr/local/devel/ANNOTATION/workflow-2.2

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
echo "To launch viewer execute... $INSTALL_DIR/run_wfmonitor.sh -i $instance"
echo "Or launch via the web http://xmen:8080/tigr-scripts/papyrus/cgi-bin/show_pipeline.cgi?xmltemplate=$instance"

exec RunWorkflow -i $instance -l $log $debugstr

