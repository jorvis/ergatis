#!/bin/sh
#This file replaces runtime vars $;ASMBL$; $;ASMBL_FILE$; $;PROJECT$; in the template .xml and .ini files



ASMBLKEY=";ASMBL;"
ASMBLFILEKEY=";ASMBL_FILE;"
SDBKEY=";SOURCE_DATABASE;"
SDBLCKEY=";SOURCE_DATABASE_LC;"
INSTANCEKEY=";INSTANCE;"
TDBKEY=";TARGET_DATABASE;"
TDBLCKEY=";TARGET_DATABASE_LC;"
TARGET_SERVER=";TARGET_SERVER;"
DATE=";DATE;"
PROJECT=";PROJECT;"
BSML_DOC=";BSML_DOCUMENT;"
BSML_PREFIX_LC=";BSML_PREFIX_LC;"
RUNPATHPREFIX=.
WORKFLOWTEMPLATEDIR=.
WRAPPERPATH=.

while getopts t:b:w:r:v:d:p:h opt
do case "$opt" in
      t) target_database=$OPTARG;;
      b) bsml_doc=$OPTARG;;
      w) wfname=$OPTARG;;
      r) wfnamedir=$OPTARG;;
      v) target_server=$OPTARG;;
      d) date=$OPTARG;;
      p) project=$OPTARG;;
      h) echo "Usage: `basename $0` -t target_database -v target_server -b bsml_doc -w wfname -r wfnamedir -d date -p project";
         echo;
	 exit 1;;
      esac
done
shift `expr $OPTIND - 1`	  

tdbkey=`echo $target_database | tr '[a-z]' '[A-Z]'`
tdblckey=`echo $target_database | tr '[A-Z]' '[a-z]'`


echo "Storing workflow xml and ini to $wfnamedir"
echo "For workflow:$wfname setting:"
echo "target_database:$target_database"
echo "bsml_doc:$bsml_doc"
echo "date:$date"
echo "target_server:$target_server"
echo "project:$project"


echo "Template xml and ini file utilized:"
echo "$WORKFLOWTEMPLATEDIR/${wfname}_template.xml"
echo "$WORKFLOWTEMPLATEDIR/$wfname.ini"

#echo "Your arguments were:"
#echo "-t $target_database"
#echo "-w $wfname"
#echo "-r $wfnamedir"
#echo "-b $bsml_doc"
#echo "-v $target_server"
#echo "-d $date"
#echo "-p $project"


cat $WORKFLOWTEMPLATEDIR/$wfname.ini | sed "s/$ASMBLKEY/$asmbl/g" | sed "s/$BSML_DOC/$bsml_doc/g" | sed "s/$PROJECT/$project/g" | sed "s/$SDBKEY/$sdbkey/g" | sed "s/$SDBLCKEY/$sdblckey/g" | sed "s/$DATE/$date/g" | sed "s/$TDBKEY/$tdbkey/g" | sed "s/$TDBLCKEY/$tdblckey/g " | sed "s/$ASMBLFILEKEY/$asmbl_file/g" | sed "s/$INSTANCEKEY/$$/g" | sed "s/$TARGET_SERVER/$target_server/g" > $wfnamedir/$wfname.ini


cat $WORKFLOWTEMPLATEDIR/${wfname}_template.xml | sed "s/$ASMBLKEY/$asmbl/g" | sed "s/$BSML_DOC/$bsml_doc/g" | sed "s/$PROJECT/$project/g" | sed "s/$SDBKEY/$sdbkey/g" | sed "s/$SDBLCKEY/$sdblckey/g" | sed "s/$DATE/$date/g" | sed "s/$TDBKEY/$tdbkey/g" | sed "s/$TDBLCKEY/$tdblckey/g" | sed "s/$ASMBLFILEKEY/$asmbl_file/g" | sed "s/$INSTANCEKEY/$$/g" | sed "s/$TARGET_SERVER/$target_server/g" > $wfnamedir/${wfname}_template.xml


