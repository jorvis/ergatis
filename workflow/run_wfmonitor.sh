#!/bin/sh

INSTALL_DIR=.

while getopts i:h opt
do case "$opt" in
      i) instance=$OPTARG;;
      h) echo "Usage: `basename $0` -i instance file"
         echo;
	 exit;;
      esac
done

source $INSTALL_DIR/exec_env.sh
WorkflowMonitor $instance  1>/dev/null 2>/dev/null
