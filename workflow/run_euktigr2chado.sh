#!/bin/sh
#-----------------------------------------------------------------------------------------------------------------------------------------
#
# run_euktigr2chado.sh
# Jay Sundaram
#
#
# This shell script will execute the euktigr2chado workflow
# 1) Accepts and verifies arguments.
# 2) source's wfenv_bash.sh (which exports workflow environmental variables)
# 3) creates workflow temporary directory
# 4) makes system call to set_runtime_vars_euktigr2chado.sh (which performs sed replacement on instance XML template and .ini files)
# 5) creates workflow instance
# 6) executes workflow instance
#
#
#-----------------------------------------------------------------------------------------------------------------------------------------
RUNPATHPREFIX=.
WORKFLOWTEMPLATEDIR=.
WRAPPERPATH=.

#------------------------------------------------------------------
# Process command options
#
# euktigr2chado expects the arguments specified below:
#
#------------------------------------------------------------------
while getopts u:p:s:t:a:v:h:f opt
do case "$opt" in
      s) source_database=$OPTARG;;
      t) target_database=$OPTARG;;
      a) asmbl_id=$OPTARG;;
      f) asmbl_id_file=$OPTARG;;
      v) target_server=$OPTARG;;
      u) username=$OPTARG;;
      p) password=$OPTARG;;
      h) echo "Usage: `basename $0` -u username -p password -s source_database -a asmbl_id|-f asmbl_id_file -t target_database -v target_server";
	  echo;
	  echo "You must specify -u, -p, -s, -a|-f, -t, and -v";
	  exit;;
  esac
done
	  

#-----------------------------------------------------------------
# Determine whether essential arguments were defined/specified
#
#-----------------------------------------------------------------

if [ -z $username ]
then
    echo "Usage: 'basename $0' -u username -p password -s source_database -a asmbl_id|-f asmbl_id_file -t target_database -v target_server";
    echo;
    echo "You must specify username";
    exit 1;
fi

if [ -z $password ]
then
    echo "Usage: 'basename $0' -u username -p password -s source_database -a asmbl_id|-f asmbl_id_file -t target_database -v target_server";
    echo;
    echo "You must specify password";
    exit 1;
fi

if [ -z $source_database ]
then
    echo "Usage: 'basename $0' -u username -p password -s source_database -a asmbl_id|-f asmbl_id_file -t target_database -v target_server";
    echo;
    echo "You must specify source_database";
    exit 1;
fi

if [ -z $asmbl_id ] && [ -z $asmbl_id_file ]
    then
    echo "Usage: 'basename $0' -u username -p password -s source_database -a asmbl_id|-f asmbl_id_file -t target_database -v target_server";
    echo;
    echo "You must specify either specify asmbl_id or asmbl_id_file";
    exit 1; 
fi

if [ -z $target_database ]
then
    echo "Usage: 'basename $0' -u username -p password -s source_database -a asmbl_id|-f asmbl_id_file -t target_database -v target_server";
    echo;
    echo "You must specify target_database";
    exit 1;
fi

if [ -z $target_server ]
then
    echo "Usage: 'basename $0' -u username -p password -s source_database -a asmbl_id|-f asmbl_id_file -t target_database -v target_server";
    echo;
    echo "You must specify target_server";
    exit 1;
fi


#--------------------------------------------------------------------------------------------
# Set up workflow vars
#
#--------------------------------------------------------------------------------------------
databasekey=`echo $source_database | tr '[a-z]' '[A-Z]'`


environ=$WRAPPERPATH/wfenv_bash.sh

#---------------------------------------------------------------------------------------------
# Check whether the environment variables exporter script exists and is executable...
#
#---------------------------------------------------------------------------------------------
if [ -x $environ ] 
    then
    source $environ
else
    echo "$environ does not exist or is not executable"
    exit 1;
fi


#---------------------------------------------------------------------------------------------
# Check whether target_server name is either SYBIL or SYBTIGR
#
#---------------------------------------------------------------------------------------------
if [  $target_server = "SYBIL" ] || [ $target_server = "SYBTIGR" ]
    then 
#    echo "$target_server is acceptable"
else
    echo "$target_server is not acceptable"
    echo "Accepted SYBASE servers are SYBIL and SYBTIGR"
    exit 1;
fi






wfname="euktigr2chado"
wfnamedir="$RUNPATHPREFIX/$databasekey/workflow/$wfname/$$"

#echo "wfname       : $wfname"
#echo "wfnamedir    : $wfnamedir"
#echo "process id   : $$"
#echo "runpathprefix: $RUNPATHPREFIX"
#echo "databasekey  : $databasekey"
#echo "path         : $PATH"
#echo "wf_root      : $WF_ROOT"

#---------------------------------------------------------------------------------------------
# Create the workflow working directory
#
#---------------------------------------------------------------------------------------------
mkdir -p $wfnamedir
if [ ! -d $wfnamedir ] 
    then
    echo "$wfnamedir was not created"
fi


#---------------------------------------------------------------------------------------------
# Check whether the asmbl_id_file exists and is readable
#
#---------------------------------------------------------------------------------------------
if [ ! -r $asmbl_id_file ] 
    then
    echo "$asmbl_id_file does not exist or does not have read permissions"
fi





#---------------------------------------------------------------------------------------------
# if asmbl_id was defined, run the variable replacement shell script with -a option
#
#---------------------------------------------------------------------------------------------
if [ -n $asmbl_id ]
    then
    $WRAPPERPATH/set_runtime_vars_euktigr2chado.sh -u $username -p $password -s $source_database -t $target_database  -w $wfname -r $wfnamedir -a $asmbl_id -v $target_server
elif [ -n $asmbl_id_file ]
    then
    #
    # if the asmbl_id_file was defined, run the variable replacement shell script with -f option
    #
    $WRAPPERPATH/set_runtime_vars_euktigr2chado.sh -u $username -p $password -s $source_database -t $target_database  -w $wfname -r $wfnamedir -f $asmbl_id_file -v $target_server
else
    echo "You must either specify the asmbl_id or the asmbl_id_file"
    exit 1;
fi

#----------------------------------------------------------------------------------------------
# Instantiate and execute the workflow instance
#
#----------------------------------------------------------------------------------------------


#/usr/local/devel/ANNOTATION/workflow/CreateWorkflow -t $wfnamedir/${wfname}_template.xml -c $wfnamedir/$wfname.ini -i $wfnamedir/$wfname.xml -l $wfnamedir/$wfname.$$.log

#/usr/local/devel/ANNOTATION/workflow/RunWorkflow -i $wfnamedir/$wfname.xml -l $wfnamedir/$wfname.$$.log

