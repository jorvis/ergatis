#!/bin/sh

KEY=";DATABASE;"
KEY2=";PROJECT;"

if [ ! -d $1 ]
then
    mkdir $1
else
   echo "The directory '$1' already exists"
   exit 1
fi

echo "Storing workflow xml and ini to $1"
echo "Database to retrieve from is $1"
echo "New project name is $2"

#echo "Saving $1/bit_score.ini"
cat db2bsml.ini | sed "s/$KEY/$1/g" | sed "s/$KEY2/$2/g" > $1/db2bsml.ini

#echo "Saving $1/bit_score_template.xml"
cat db2bsml_template.xml | sed "s/$KEY/$1/g" > $1/db2bsml_template.xml

echo "Start workflow with cd $1;"
echo "RunWorkflow -t db2bsml_template.xml -c db2bsml.ini -i db2bsml.$1.xml -l db2bsml.$1.log"
