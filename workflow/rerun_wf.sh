#!/bin/sh

while getopts i:l:h opt
do case "$opt" in
      i) instance=$OPTARG;;
      l) log=$OPTARG;;
      h) echo "Usage: `basename $0` -i instance file -l log file"
         echo;
	 exit;;
      esac
done

source /export/CVS/papyrus/workflow/distrib-promer/exec_env.sh

echo "Running workflow instance $instance"
echo "To launch viewer execute... /export/CVS/papyrus/workflow/distrib-promer/run_wfmonitor.sh -i $instance"
RunWorkflow -i $instance -l $log 1>>shell_log.txt 2>>shell_log.txt
