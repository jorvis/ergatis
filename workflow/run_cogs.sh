#!/bin/sh

#####################################################
#
# run_cogs.sh - interface for the cogs analysis workflow
#
# Args:
#   -d Database Name - to direct to BSML Repository
#   -s search [all_vs_all or blastp]
#
# ./run_cogs.sh -d TRYP -s all_vs_all

RUNPATH=.

#process command options

while getopts d:s:h opt
do case "$opt" in
      d) database=$OPTARG;;
      s) search=$OPTARG;;
      h) echo "Usage: `basename $0` -d dbname -s search";
	 echo "where search = [all_vs_all | blastp]";
	  exit;;
      esac
done

#Set up workflow vars
source $RUNPATH/wfenv_bash.sh

project="$database$$"
program="cogs"

    $RUNPATH/set_runtime_vars.sh -d $database -p $project -r $program -a $search

/usr/local/devel/ANNOTATION/workflow/CreateWorkflow -t $project/${program}_template.xml -c $project/$program.ini -i $project/$program.xml -l $project/$program.$$.log
/usr/local/devel/ANNOTATION/workflow/RunWorkflow -i $project/$program.xml -l $project/$program.$$.log
