#!/bin/sh


RUNPATHPREFIX=.
WORKFLOWTEMPLATEDIR=.
WRAPPERPATH=.

#process command options

while getopts t:b:v:d:p:h opt
  do case "$opt" in
      t) target_database=$OPTARG;;
      b) bsml_doc=$OPTARG;;
      v) target_server=$OPTARG;;
      d) date=$OPTARG;;
      p) project=$OPTARG;;
      h) echo "Usage: `basename $0` -b bsml_doc -t target_database -v target_server -d date -p project";
	  echo;
	  exit;;
  esac
done
	  

if [ -z $target_database ]
then
    echo "Usage: 'basename $0' -b bsml_doc -t target_database -v target_server -d date -p project";
    echo;
    exit 1;
fi

if [ -z $bsml_doc ]
then
    echo "Usage: 'basename $0' -b bsml_doc -t target_database -v target_server -d date -p project";
    echo;
    exit 1;
fi

if [ -z $target_server ]
then
    echo "Usage: 'basename $0' -b bsml_doc -t target_database -v target_server -d date -p project";
    echo;
    exit 1;
fi

if [ -z $date ]
then
    echo "Usage: 'basename $0' -b bsml_doc -t target_database -v target_server -d date -p project";
    echo;
    exit 1;
fi


#Set up workflow vars
databasekey=`echo $database | tr '[a-z]' '[A-Z]'`
source $WRAPPERPATH/wfenv_bash.sh

wfname="bsmlCOG2Chado"
wfnamedir="$RUNPATHPREFIX/$databasekey/workflow/$wfname/$$"

#echo "wfname: $wfname"
#echo "wfnamedir: $wfnamedir"
#echo "process id: $$"

mkdir -p $wfnamedir
$WRAPPERPATH/set_runtime_vars_bsmlCOG2Chado.sh -t $target_database -w $wfname -r $wfnamedir -b $bsml_doc -v $target_server -d $date -p $project

#/usr/local/devel/ANNOTATION/workflow/CreateWorkflow -t $wfnamedir/${wfname}_template.xml -c $wfnamedir/$wfname.ini -i $wfnamedir/$wfname.xml -l $wfnamedir/$wfname.$$.log
#/usr/local/devel/ANNOTATION/workflow/RunWorkflow -i $wfnamedir/$wfname.xml -l $wfnamedir/$wfname.$$.log

