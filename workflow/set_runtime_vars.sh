#!/bin/sh
#This file replaces runtime vars $;ASMBL$; $;ASMBL_FILE$; $;PROJECT$;



ASMBLKEY=";ASMBL;"
ASMBLFILEKEY=";ASMBL_FILE;"
DBKEY=";DATABASE;"
DBLCKEY=";DATABASE_LC;"
PROJECTKEY=";PROJECT;"

while getopts d:a:f:p:r:h opt
do case "$opt" in
      d) database=$OPTARG;;
      a) asmbl=$OPTARG;;
      f) asmbl_file=$OPTARG;;
      p) project=$OPTARG;;
      r) program=$OPTARG;;
      h) echo "Usage: `basename $0` -d dbname -p project_name -r program [-a asmbl_id] [-f asmbl_file]";
          echo;
	  echo "You must specify -a or -f";
	  exit 1;;
      esac
done
shift `expr $OPTIND - 1`	  

if [ ! -d $project ]
then
    mkdir $project
fi

databasekey=`echo $database | tr '[a-z]' '[A-Z]'`
databaselckey=`echo $database | tr '[A-Z]' '[a-z]'`

echo "Storing workflow xml and ini to $project"
echo "Setting asmbl:$asmbl asmbl_file:$asmbl_file database:$database for $program"
cat $program.ini | sed "s/$PROJECTKEY/$project/g" | sed "s/$ASMBLKEY/$asmbl/g" | sed "s/$DBKEY/$databasekey/g" | sed "s/$DBLCKEY/$databaselckey/g" | sed "s/$ASMBLFILEKEY/$asmbl_file/g" > $project/$program.ini
cat ${program}_template.xml | sed "s/$ASMBLKEY/$asmbl/g" | sed "s/$DBKEY/$databasekey/g" | sed "s/$DBLCKEY/$databaselckey/g" | sed "s/$ASMBLFILEKEY/$asmbl_file/g" > $project/${program}_template.xml



