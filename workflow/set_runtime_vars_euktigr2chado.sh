#!/bin/sh
#----------------------------------------------------------------------------------------------------
#
#
#
#
#
#
#
# This file replaces runtime vars $;ASMBL$; $;ASMBL_FILE$; $;PROJECT$; in the template xml and 
#
#
#----------------------------------------------------------------------------------------------------

ASMBLKEY=";ASMBL;"
ASMBLFILEKEY=";ASMBL_FILE;"
SDBKEY=";SOURCE_DATABASE;"
SDBLCKEY=";SOURCE_DATABASE_LC;"
INSTANCEKEY=";INSTANCE;"
TDBKEY=";TARGET_DATABASE;"
TDBLCKEY=";TARGET_DATABASE_LC;"
TARGET_SERVER=";TARGET_SERVER;"
USERNAME=";USERNAME;"
PASSWORD=";PASSWORD;"

RUNPATHPREFIX=.
WORKFLOWTEMPLATEDIR=.
WRAPPERPATH=.

#---------------------------------------------------------------------------------
# Accept and verify arguments
#
#---------------------------------------------------------------------------------
while getopts u:p:s:t:a:f:w:r:v:h opt
  do case "$opt" in
  u) username=$OPTARG;;
      p) password=$OPTARG;;
      s) source_database=$OPTARG;;
      t) target_database=$OPTARG;;
      v) target_server=$OPTARG;;
      a) asmbl_id=$OPTARG;;
      f) asmbl_id_file=$OPTARG;;
      w) wfname=$OPTARG;;
      r) wfnamedir=$OPTARG;;
      h) echo "Usage: `basename $0` -u username -p password -s source_database -t target_database -a asmbl_id|-f asmbl_id_file -t target_database -v target_server -w wfname -r wfnamedir";
	  echo;
	  exit 1;;
  esac
done
shift `expr $OPTIND - 1`	  

if [ -z $username ]
    then
    echo "Usage: 'basename $0' -u username -p password -s source_database -t target_database -a asmbl_id|-f asmbl_id_file -t target_database -v target_server -w wfname -r wfnamedir";
    echo;
    echo "You must specify -u username";
    exit 1;
fi

if [ -z $password ]
    then
    echo "Usage: 'basename $0' -u username -p password -s source_database -t target_database -a asmbl_id|-f asmbl_id_file -t target_database -v target_server -w wfname -r wfnamedir";
    echo;
    echo "You must specify -p password";
    exit 1;
fi

if [ -z $source_database ]
    then
    echo "Usage: 'basename $0' -u username -p password -s source_database -t target_database -a asmbl_id|-f asmbl_id_file -t target_database -v target_server -w wfname -r wfnamedir";
    echo;
    echo "You must specify -s source_database";
    exit 1;
fi

if [ -z $target_database ]
    then
    echo "Usage: 'basename $0' -u username -p password -s source_database -t target_database -a asmbl_id|-f asmbl_id_file -t target_database -v target_server -w wfname -r wfnamedir";
    echo;
    echo "You must specify -t target_database";
    exit 1;
fi

if [ -z $asmbl_id ] &&  [ -z $asmbl_id_file ]
    then
    echo "Usage: 'basename $0' -u username -p password -s source_database -t target_database -a asmbl_id|-f asmbl_id_file -t target_database -v target_server -w wfname -r wfnamedir";
    echo;
    echo "You must specify either asmbl_id or asmbl_id_file";
    exit 1;
fi

if [ -z $target_server ]
    then
    echo "Usage: 'basename $0' -u username -p password -s source_database -t target_database -a asmbl_id|-f asmbl_id_file -t target_database -v target_server -w wfname -r wfnamedir";
    echo;
    echo "You must specify -v target_server";
    exit 1;
fi

if [ -z $wfname ]
    then
    echo "Usage: 'basename $0' -u username -p password -s source_database -t target_database -a asmbl_id|-f asmbl_id_file -t target_database -v target_server -w wfname -r wfnamedir";
    echo;
    echo "You must specify -w wfname";
    exit 1;
fi

if [ -z $wfnamedir ]
    then
    echo "Usage: 'basename $0' -u username -p password -s source_database -t target_database -a asmbl_id|-f asmbl_id_file -t target_database -v target_server -w wfname -r wfnamedir";
    echo;
    echo "You must specify -r wfnamedir";
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

#---------------------------------------------------------------------------------------------
# Check whether wfnamedir exists
#
#---------------------------------------------------------------------------------------------
if [  ! -d $wfnamedir ]
    then 
    echo "$wfnamedir does not exist"
    exit 1;
fi


sdbkey=`echo $source_database | tr '[a-z]' '[A-Z]'`
sdblckey=`echo $source_database | tr '[A-Z]' '[a-z]'`
tdbkey=`echo $target_database | tr '[a-z]' '[A-Z]'`
tdblckey=`echo $target_database | tr '[A-Z]' '[a-z]'`



#echo "Setting following variables:"
#echo "source legacy database:$source_database"
#echo "target chado database :$target_database"
#echo "asmbl_id              :$asmbl_id"
#echo "wfname                :$wfname"
#echo "wfnamedir             :$wfnamedir"
#echo "" 
echo "Storing workflow xml and ini to $wfnamedir"
echo "Template xml and ini file utilized:"
echo "$WORKFLOWTEMPLATEDIR/${wfname}_template.xml"
echo "$WORKFLOWTEMPLATEDIR/$wfname.ini"
echo ""
#echo "Your arguments were:"
#echo "-u $username"
#echo "-p $password"
#echo "-s $source_database"
#echo "-t $target_database"
#echo "-v $target_server"
#echo "-w $wfname"
#echo "-r $wfnamedir"
#echo "-a $asmbl_id"
#echo "-f $asmbl_id_file"


cat $WORKFLOWTEMPLATEDIR/$wfname.ini | sed "s/$USERNAME/$username/g" | sed "s/$PASSWORD/$password/g" | sed "s/$ASMBLKEY/$asmbl_id/g" | sed "s/$SDBKEY/$sdbkey/g" | sed "s/$SDBLCKEY/$sdblckey/g" | sed "s/$TDBKEY/$tdbkey/g" | sed "s/$TDBLCKEY/$tdblckey/g " | sed "s/$ASMBLFILEKEY/$asmbl_id_file/g" | sed "s/$INSTANCEKEY/$$/g" | sed "s/$TARGET_SERVER/$target_server/g" > $wfnamedir/$wfname.ini
cat $WORKFLOWTEMPLATEDIR/${wfname}_template.xml | sed "s/$USERNAME/$username/g" | sed "s/$PASSWORD/$password/g" | sed "s/$ASMBLKEY/$asmbl_id/g" | sed "s/$SDBKEY/$sdbkey/g" | sed "s/$SDBLCKEY/$sdblckey/g" | sed "s/$TDBKEY/$tdbkey/g" | sed "s/$TDBLCKEY/$tdblckey/g" | sed "s/$ASMBLFILEKEY/$asmbl_id_file/g" | sed "s/$INSTANCEKEY/$$/g" | sed "s/$TARGET_SERVER/$target_server/g" > $wfnamedir/${wfname}_template.xml


