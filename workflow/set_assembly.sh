#!/bin/sh

KEY=";ASMBL;"

if [ ! -d $1 ]
then
    mkdir $1
fi

echo "Storing workflow xml and ini to $1"

echo "Saving $1/bit_score.ini"
cat bit_score.ini | sed "s/$KEY/$1/" > $1/bit_score.ini
echo "Saving $1/pe.ini"
cat pe.ini | sed "s/$KEY/$1/" > $1/pe.ini
echo "Saving $1/all_vs_all.ini"
cat all_vs_all.ini | sed "s/$KEY/$1/" > $1/all_vs_all.ini

echo "Saving $1/bit_score_template.xml"
cat bit_score_template.xml | sed "s/$KEY/$1/" > $1/bit_score_template.xml
echo "Saving $1/pe_template.xml"
cat pe_template.xml | sed "s/$KEY/$1/" > $1/pe_template.xml
echo "Saving $1/all_vs_all_template.xml"
cat all_vs_all_template.xml | sed "s/$KEY/$1/" > $1/all_vs_all_template.xml


echo "Start workflow with cd $1;"
echo "RunWorkflow -t all_vs_all_template.xml -c all_vs_all.ini -i allvsall.$1.xml -l allvsall.$1.log"
