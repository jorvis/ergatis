#!/bin/sh

KEY=";ASMBL;"
DB=";PROJECT;"



if [ ! -d $1 ]
then
    mkdir $1
fi

echo "Storing workflow xml and ini to $1"
echo "new asmbl_id is $1"
echo "new project name is $2"



echo "Saving $1/bit_score.ini"
cat bit_score.ini | sed "s/$KEY/$1/g" | sed "s/$DB/$2/g" > $1/bit_score.ini

echo "Saving $1/pe.ini"
cat pe.ini | sed "s/$KEY/$1/g" | sed "s/$DB/$2/g" > $1/pe.ini

echo "Saving $1/all_vs_all.ini"
cat all_vs_all.ini | sed "s/$KEY/$1/g" | sed "s/$DB/$2/g" > $1/all_vs_all.ini

echo "Saving $1/all_vs_all_fast.ini"
cat all_vs_all_fast.ini | sed "s/$KEY/$1/g" | sed "s/$DB/$2/g" > $1/all_vs_all_fast.ini

echo "Saving $1/moaf.ini"
cat moaf.ini | sed "s/$KEY/$1/g" | sed "s/$DB/$2/g" > $1/moaf.ini

echo "Saving $1/bit_score_template.xml"
cat bit_score_template.xml | sed "s/$KEY/$1/g" > $1/bit_score_template.xml
echo "Saving $1/pe_template.xml"
cat pe_template.xml | sed "s/$KEY/$1/g" > $1/pe_template.xml
echo "Saving $1/all_vs_all_template.xml"
cat all_vs_all_template.xml | sed "s/$KEY/$1/g" > $1/all_vs_all_template.xml
echo "Saving $1/all_vs_all_fast_template.xml"
cat all_vs_all_fast_template.xml | sed "s/$KEY/$1/g" > $1/all_vs_all_fast_template.xml

echo "Saving $1/moaf_template.xml"
cat moaf_template.xml | sed "s/$KEY/$1/g" > $1/moaf_template.xml


echo "Start workflow with cd $1;"
echo "RunWorkflow -t moaf_template.xml -c moaf.ini -i moaf.$1.xml -l moaf.$1.log"
echo "RunWorkflow -t all_vs_all_template.xml -c all_vs_all.ini -i allvsall.$1.xml -l allvsall.$1.log"
echo "RunWorkflow -t all_vs_all_fast_template.xml -c all_vs_all_fast.ini -i allvsall_fast.$1.xml -l allvsall_fast.$1.log"
echo "RunWorkflow -t pe_template.xml -c pe.ini -i pe.$1.xml -l pe.$1.log"
echo "RunWorkflow -t bit_score_template.xml -c bit_score.ini -i bit_score.$1.xml -l bit_score.$1.log"
