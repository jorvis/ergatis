#!/bin/sh


RUNPATHPREFIX=.
WORKFLOWTEMPLATEDIR=.
WRAPPERPATH=.

#process command options

while getopts s:t:v:h opt
do case "$opt" in
      s) source=$OPTARG;;
      t) target=$OPTARG;;
      v) target_server=$OPTARG;;
      h) echo "Usage: `basename $0` -s source_database_name -t target_database_name -v target_server";
          echo;
	  echo "You must specify -s, -t and -v";
	  exit;;
      esac
done
	  
if [ -z $source ]
then
    echo "Usage: 'basename $0' -s source_database_name dbname -t target_database_name -v target_server";
    echo;
    echo "You must specify -s, -t and -v";
    exit 1;
fi

if [ -z $target ]
then
    echo "Usage: 'basename $0' -s source_database_name dbname -t target_database_name -v target_server";
    echo;
    echo "You must specify -s, -t and -v";
    exit 1;
fi

if [ -z $target_server ]
then
    echo "Usage: 'basename $0' -s source_database_name dbname -t target_database_name -v target_server";
    echo;
    echo "You must specify -s, -t and -v";
    exit 1;
fi



#Set up workflow vars
databasekey=`echo $database | tr '[a-z]' '[A-Z]'`
source $WRAPPERPATH/wfenv_bash.sh

wfname="proktigr2chado"
wfnamedir="$RUNPATHPREFIX/$databasekey/workflow/$wfname/$$"

#echo "wfname: $wfname"
#echo "wfnamedir: $wfnamedir"
#echo "process id: $$"


mkdir -p $wfnamedir
$WRAPPERPATH/set_runtime_proktigr2chado_vars.sh -s $source -t $target -v $target_server -w $wfname -r $wfnamedir

#/usr/local/devel/ANNOTATION/workflow/CreateWorkflow -t $wfnamedir/${wfname}_template.xml -c $wfnamedir/$wfname.ini -i $wfnamedir/$wfname.xml -l $wfnamedir/$wfname.$$.log
#/usr/local/devel/ANNOTATION/workflow/RunWorkflow -i $wfnamedir/$wfname.xml -l $wfnamedir/$wfname.$$.log

