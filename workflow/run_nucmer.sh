#!/bin/sh

RUNPATHPREFIX=.
WORKFLOWTEMPLATEDIR=.
WRAPPERPATH=.

#process command options

while getopts d:a:f:h opt
do case "$opt" in
      d) database=$OPTARG;;
      a) asmbl=$OPTARG;;
      f) asmbl_file=$OPTARG;;
      h) echo "Usage: `basename $0` -d dbname [-a asmbl_id] [-f asmbl_file]";
          echo;
	  echo "You must specify -a or -f";
	  exit;;
      esac
done
	  
if [ -z $database ]
then
    echo "Usage: 'basename $0' -d dbname [-a asmbl_id] [-f asmbl_file]";
    echo;
    echo "You must specify -a or -f";
    exit 1;
fi

#Set up workflow vars
databasekey=`echo $database | tr '[a-z]' '[A-Z]'`
source $WRAPPERPATH/wfenv_bash.sh

if [ "$asmbl_file" ]
then
    wfname="nucmer_multi"
    wfnamedir="$RUNPATHPREFIX/$databasekey/workflow/$wfname/$$"
    mkdir -p $wfnamedir
    $WRAPPERPATH/set_runtime_vars.sh -f $asmbl_file -d $database -w $wfname -r $wfnamedir -a $asmbl
    cp $asmbl_file $wfnamedir
elif [ "$asmbl" ]
then 
    wfname="nucmer"
    wfnamedir="$RUNPATHPREFIX/$databasekey/workflow/$wfname/$$"
    mkdir -p $wfnamedir
    $WRAPPERPATH/set_runtime_vars.sh -a $asmbl -d $database  -w $wfname -r $wfnamedir 
else
    echo "You must specify -a or -f"
    exit;
fi
/usr/local/devel/ANNOTATION/workflow/CreateWorkflow -t $wfnamedir/${wfname}_template.xml -c $wfnamedir/$wfname.ini -i $wfnamedir/$wfname.xml -l $wfnamedir/$wfname.$$.log
/usr/local/devel/ANNOTATION/workflow/RunWorkflow -i $wfnamedir/$wfname.xml -l $wfnamedir/$wfname.$$.log
