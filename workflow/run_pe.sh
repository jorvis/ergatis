#!/bin/sh

RUNPATH=.

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
source $RUNPATH/wfenv_bash.sh

project="$database$$"

if [ "$asmbl_file" ]
then
    program="pe_multi"
    $RUNPATH/set_runtime_vars.sh -f $asmbl_file -d $database -p $project -r $program -a $asmbl
    cp $asmbl_file $project
elif [ "$asmbl" ]
then 
    program="pe"
     $RUNPATH/set_runtime_vars.sh -a $asmbl -d $database -p $project -r $program
else
    echo "You must specify -a or -f"
    exit;
fi
/usr/local/devel/ANNOTATION/workflow/CreateWorkflow -t $project/${program}_template.xml -c $project/$program.ini -i $project/$program.xml -l $project/$program.$$.log
/usr/local/devel/ANNOTATION/workflow/RunWorkflow -i $project/$program.xml -l $project/$program.$$.log
