#!/bin/sh
#This file replaces runtime vars $;ASMBL$; $;ASMBL_FILE$; $;PROJECT$; in the template xml and 



ASMBLKEY=";ASMBL;"
ASMBLFILEKEY=";ASMBL_FILE;"
SDBKEY=";SOURCE_DATABASE;"
SDBLCKEY=";SOURCE_DATABASE_LC;"
INSTANCEKEY=";INSTANCE;"
TDBKEY=";TARGET_DATABASE;"
TDBLCKEY=";TARGET_DATABASE_LC;"
TARGET_SERVER=";TARGET_SERVER;"
RUNPATHPREFIX=.
WORKFLOWTEMPLATEDIR=.
WRAPPERPATH=.

while getopts s:t:a:f:w:r:v:h opt
do case "$opt" in
      s) database=$OPTARG;;
      t) target_database=$OPTARG;;
      a) asmbl=$OPTARG;;
      f) asmbl_file=$OPTARG;;
      w) wfname=$OPTARG;;
      r) wfnamedir=$OPTARG;;
      v) target_server=$OPTARG;;
      h) echo "Usage: `basename $0` -s database -t target_database -v target_server -w wf_name -r run_dir";
         echo;
	 exit 1;;
      esac
done
shift `expr $OPTIND - 1`	  

sdbkey=`echo $database | tr '[a-z]' '[A-Z]'`
sdblckey=`echo $database | tr '[A-Z]' '[a-z]'`
tdbkey=`echo $target_database | tr '[a-z]' '[A-Z]'`
tdblckey=`echo $target_database | tr '[A-Z]' '[a-z]'`



echo "Storing workflow xml and ini to $wfnamedir"
echo "Setting database:$database target_database:$target_database for $wfname"

echo "Template xml and ini file utilized:"
echo "$WORKFLOWTEMPLATEDIR/${wfname}_template.xml"
echo "$WORKFLOWTEMPLATEDIR/$wfname.ini"

#echo "Your arguments were:"
#echo "-s $database"
#echo "-t $target_database"
#echo "-w $wfname"
#echo "-r $wfnamedir"
#echo "-v $target_server"

cat $WORKFLOWTEMPLATEDIR/$wfname.ini | sed "s/$ASMBLKEY/$asmbl/g" | sed "s/$SDBKEY/$sdbkey/g" | sed "s/$SDBLCKEY/$sdblckey/g" | sed "s/$TDBKEY/$tdbkey/g" | sed "s/$TDBLCKEY/$tdblckey/g " | sed "s/$ASMBLFILEKEY/$asmbl_file/g" | sed "s/$INSTANCEKEY/$$/g" | sed "s/$TARGET_SERVER/$target_server/g"  > $wfnamedir/$wfname.ini
cat $WORKFLOWTEMPLATEDIR/${wfname}_template.xml | sed "s/$ASMBLKEY/$asmbl/g" | sed "s/$SDBKEY/$sdbkey/g" | sed "s/$SDBLCKEY/$sdblckey/g" | sed "s/$TDBKEY/$tdbkey/g" | sed "s/$TDBLCKEY/$tdblckey/g" | sed "s/$ASMBLFILEKEY/$asmbl_file/g" | sed "s/$INSTANCEKEY/$$/g" | sed "s/$TARGET_SERVER/$target_server/g" > $wfnamedir/${wfname}_template.xml


