#!/bin/sh

#####################################################
#
# run_snps.sh - interface for the snps analysis workflow
#
# Args:
#   -d Database Name - to direct to BSML Repository
#   -a Reference Assembly ID
#   -f File Containing the reference id (this should be changed when ready for production)
#   -g File Containing the query ids
# ./run_snps.sh -d TB -a gmt_3810_assembly -q query_list.txt -r reference_list.txt
# 
# USAGE NOTE - the SED replacements used in set_runtime_vars.sh require special characters to be escaped...
# This needs to be addressed.
#
# ./run_snps.sh -d CHADO_TEST -q '\/usr\/local\/annotation\/CHADO_TEST\/BSML_repository\/query.txt' -r '\/usr\/local\/annotation\/CHADO_TEST\/BSML_repository\/reference.txt' 
# 

RUNPATHPREFIX=.
WORKFLOWTEMPLATEDIR=.
WRAPPERPATH=.

#process command options

while getopts d:a:q:r:h opt
do case "$opt" in
      d) database=$OPTARG;;
      a) refid=$OPTARG;;
      q) query=$OPTARG;;
      r) ref=$OPTARG;;
      h) echo "Usage: `basename $0` -d dbname -a reference assembly id -f reference assembly list -g query assembly list";
	  exit;;
      esac
done

#Set up workflow vars
databasekey=`echo $database | tr '[a-z]' '[A-Z]'`
source $WRAPPERPATH/wfenv_bash.sh

wfname="snp"
wfnamedir="$RUNPATHPREFIX/$databasekey/workflow/$wfname/$$"
mkdir -p $wfnamedir

$WRAPPERPATH/set_runtime_vars.sh -d $database -w $wfname -r $wfnamedir -f $ref -g $query -a $refid

/usr/local/devel/ANNOTATION/workflow/CreateWorkflow -t $wfnamedir/${wfname}_template.xml -c $wfnamedir/$wfname.ini -i $wfnamedir/$wfname.xml -l $wfnamedir/$wfname.$$.log
/usr/local/devel/ANNOTATION/workflow/RunWorkflow -i $wfnamedir/$wfname.xml -l $wfnamedir/$wfname.$$.log
