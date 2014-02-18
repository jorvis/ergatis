#!/bin/sh
# USAGE : sh download_tbl2asn.sh <path to output directory>

if [ $# -eq 1 ] && [ -d $1 ] && [ -e $1 ];
then
	OUTDIR=$1
	FILE=linux.tbl2asn.gz

	wget -t 50 -P $OUTDIR ftp://ftp.ncbi.nih.gov/toolbox/ncbi_tools/converters/by_program/tbl2asn/linux.tbl2asn.gz 2>$OUTDIR/tbl2asn.log

	if [ -s $OUTDIR/$FILE ] && [ $? -eq 0 ]; 
	then
		echo "INFO : cd $OUTDIR" >> $OUTDIR/tbl2asn.log
		cd $OUTDIR
		if [ -s tbl2asn ];
		then
			echo "INFO : rm tbl2asn if exists tbl2asn" >> tbl2asn.log
			rm tbl2asn
		fi
		echo "INFO : gunzip -d $FILE" >> tbl2asn.log
		gunzip -d $FILE
		echo "INFO : chmod +x linux.tbl2asn" >> tbl2asn.log
		chmod +x linux.tbl2asn
		echo "INFO : mv linux.tbl2asn tbl2asn" >> tbl2asn.log
		mv linux.tbl2asn tbl2asn
	else
		echo "ERROR : Unable to obtain linux.tbl2asn.gz from NCBI FTP site." >> $OUTDIR/tbl2asn.log
	fi

	if [ -s tbl2asn ];
	then
		echo "INFO : Latest version of tbl2asn is available in $OUTDIR" >> tbl2asn.log
	else 
		echo "ERROR : Failed to get latest version of tbl2asn in $OUTDIR. Check log file" >> $OUTDIR/tbl2asn.log
	fi
else 
	echo "USAGE : sh $0 <path to output direcroty>"	
fi	
