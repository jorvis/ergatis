#!/bin/sh
#This file replaces runtime vars $;ASMBL$; $;ASMBL_FILE$; $;PROJECT$; in the template xml and 



ASMBLKEY=";ASMBL;"
ASMBLFILEKEY=";ASMBL_FILE;"
DBKEY=";DATABASE;"
DBLCKEY=";DATABASE_LC;"
INSTANCEKEY=";INSTANCE;"
RUNPATHPREFIX=.
WORKFLOWTEMPLATEDIR=.
WRAPPERPATH=.

while getopts d:a:f:w:r:h opt
do case "$opt" in
      d) database=$OPTARG;;
      a) asmbl=$OPTARG;;
      f) asmbl_file=$OPTARG;;
      w) wfname=$OPTARG;;
      r) wfnamedir=$OPTARG;;
      h) echo "Usage: `basename $0` -d dbname -w wf_name -r run_dir [-a asmbl_id] [-f asmbl_file]";
          echo;
	  echo "You must specify -a or -f";
	  exit 1;;
      esac
done
shift `expr $OPTIND - 1`	  

databasekey=`echo $database | tr '[a-z]' '[A-Z]'`
databaselckey=`echo $database | tr '[A-Z]' '[a-z]'`

echo "Storing workflow xml and ini to $wfnamedir"
echo "Setting asmbl:$asmbl asmbl_file:$asmbl_file database:$database for $wfname"
cat $WORKFLOWTEMPLATEDIR/$wfname.ini | sed "s/$ASMBLKEY/$asmbl/g" | sed "s/$DBKEY/$databasekey/g" | sed "s/$DBLCKEY/$databaselckey/g" | sed "s/$ASMBLFILEKEY/$asmbl_file/g" | sed "s/$INSTANCEKEY/$$/g" > $wfnamedir/$wfname.ini
cat $WORKFLOWTEMPLATEDIR/${wfname}_template.xml | sed "s/$ASMBLKEY/$asmbl/g" | sed "s/$DBKEY/$databasekey/g" | sed "s/$DBLCKEY/$databaselckey/g" | sed "s/$ASMBLFILEKEY/$asmbl_file/g" | sed "s/$INSTANCEKEY/$$/g" > $wfnamedir/${wfname}_template.xml


