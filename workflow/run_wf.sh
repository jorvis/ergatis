#!/bin/sh

INSTALL_DIR=.

while getopts d:t:c:i:l:h opt
do case "$opt" in
      d) directory=$OPTARG;;
      t) template=$OPTARG;;
      c) ini=$OPTARG;;
      i) instance=$OPTARG;;
      l) log=$OPTARG;;
      h) echo "Usage: `basename $0` -d run_directory -t template file -c ini file -i instance file -l log file"
         echo;
	 exit;;
      esac
done
	  

source $INSTALL_DIR/exec_env.sh
cd $directory

echo "Creating workflow instance $instance"
CreateWorkflow -t $template -c $ini -i $instance -l $log 1>shell_log.txt 2>shell_log.txt

echo "Running workflow instance $instance"
echo "To launch viewer execute... $INSTALL_DIR/run_wfmonitor.sh -i $instance"

#RunWorkflow -i $instance -l $log 1>>shell_log.txt 2>>shell_log.txt &

