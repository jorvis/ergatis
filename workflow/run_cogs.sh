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

RUNPATHPREFIX=.
WORKFLOWTEMPLATEDIR=.
WRAPPERPATH=.

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

if [ -z $database ]
then
    echo "Usage: `basename $0` -d dbname -s search";
    echo "where search = [all_vs_all | blastp]";
    echo;
    echo "You must specify -d and -s";
    exit 1;
fi


#Set up workflow vars
#Set up workflow vars
databasekey=`echo $database | tr '[a-z]' '[A-Z]'`
source $WRAPPERPATH/wfenv_bash.sh

wfname="cogs"
wfnamedir="$RUNPATHPREFIX/$databasekey/workflow/$wfname/$$"
mkdir -p $wfnamedir

$WRAPPERPATH/set_runtime_vars.sh -d $database -w $wfname -r $wfnamedir -a $search

/usr/local/devel/ANNOTATION/workflow/CreateWorkflow -t $wfnamedir/${wfname}_template.xml -c $wfnamedir/$wfname.ini -i $wfnamedir/$wfname.xml -l $wfnamedir/$wfname.$$.log
/usr/local/devel/ANNOTATION/workflow/RunWorkflow -i $wfnamedir/$wfname.xml -l $wfnamedir/$wfname.$$.log

