# This program aligns the reference genome with the sequence provided.
# The indexing should be completed before this program is executed.

#!/bin/sh
# AWK has been used to parse the file.
# Author: Ashwinkumar Ganesan.
# Date Created: 01-03-2011

# This is for testing purposes.
# set -vx

# Following is the list of parameters to the shell script.
# 1  - Is the File path provided by the user or a given pipeline and the reads file in FASTQ format
# 2  - Is the location of the bwa executable file.
# 3  - Is the Output directory location.
# 4  - Mismatch Penalty.
# 5  - Maximum Open Gaps.
# 6  - Maximum Gap Extensions.
# 7  - Gap Open Penalty.
# 8  - Gap Extension Penalty.
# 9  - Maximum Number Of Threads.
# 10 - This is for the option to use BWASW.
# 11 - Number of alignments for each read.


# Copy the reference file index to the temporary location. This is because all the jobs will try to access the same file.
#cp $1* $3

# Reference File Location.
# This is the whole location.
REF_FILENAME=`echo "$1" | awk -F "|" '{print $2}'`
echo $REF_FILENAME
REF_FILENAME_WITHOUT_LOC=`echo "$REF_FILENAME" | awk -F"/" '{print $NF}'`

# First Short reads File Location.
COM_FILENAME=`echo "$1" | awk -F "|" '{print $1}'`
echo $COM_FILENAME

# This is to find the path only.
COM_FILENAME_LEN=`expr length $COM_FILENAME`
COM_FILENAME_WITHOUT_LOC=`echo $COM_FILENAME | awk -F"/" '{print $NF}'`
COM_FILENAME_WITHOUT_LOC_LEN=`expr length $COM_FILENAME_WITHOUT_LOC`
END_LOC=`expr $COM_FILENAME_LEN - $COM_FILENAME_WITHOUT_LOC_LEN`
COM_FILE_PATH=`echo $COM_FILENAME | cut -c1-${END_LOC}`

IF_PAIRED_OR_NOT=`ls -l ${COM_FILE_PATH}/*_1.fastq | wc -l`
echo ${10}
echo $IF_PAIRED_OR_NOT
echo $REF_FILENAME
echo $COM_FILENAME
echo $COM_FILE_PATH
echo $COM_FILENAME_WITHOUT_LOC

# This is a temporary piece of code that can be removed later.
if [ $IF_PAIRED_OR_NOT -eq 0 ]
then
	exit
fi

# Once the index file is created, they need to be aligned with the short 
# reads file.
# The .sai file is stored in the same output directory.
if [ ${10} -eq 0 ]
then
	if [ $IF_PAIRED_OR_NOT -eq 0 ]
	then
		$2/bwa aln $REF_FILENAME -M $4 -o $5 -e $6 -O $7 -E ${8} -t ${9} "${COM_FILE_PATH}/${COM_FILENAME_WITHOUT_LOC}.fastq" > "$3/${REF_FILENAME_WITHOUT_LOC}_${COM_FILENAME_WITHOUT_LOC}_aln_sa.sai"
		echo "File Alignment is done"

	# In case there is a paired end read, the alignment needs to be done twice.
	# Second alignment is with the seocnd short read file provided.
	elif [ $IF_PAIRED_OR_NOT -eq 1 ]
	then
   		$2/bwa aln $REF_FILENAME -M $4 -o $5 -e $6 -O $7 -E ${8} -t ${9} "${COM_FILE_PATH}/${COM_FILENAME_WITHOUT_LOC}_1.fastq" > "$3/${REF_FILENAME_WITHOUT_LOC}_${COM_FILENAME_WITHOUT_LOC}_aln_sa1.sai"
		echo "File Alignment of First File Completed"
		$2/bwa aln $REF_FILENAME -M $4 -o $5 -e $6 -O $7 -E ${8} -t ${9} "${COM_FILE_PATH}/${COM_FILENAME_WITHOUT_LOC}_2.fastq" > "$3/${REF_FILENAME_WITHOUT_LOC}_${COM_FILENAME_WITHOUT_LOC}_aln_sa2.sai"
		echo "File Alignment of Second File Completed"
	fi

	# To create the SAM file.
	# The $11 value checks if the flgag for single or paired-end read is set.
	# Accordingly the same file created.

	# "samse" option is for single reads.
	# "sampe" option is for paired reads.

	if [ $IF_PAIRED_OR_NOT -eq 0 ]
	then
		$2/bwa samse -n ${11} $REF_FILENAME "$3/${REF_FILENAME_WITHOUT_LOC}_${COM_FILENAME_WITHOUT_LOC}_aln_sa.sai" "${COM_FILE_PATH}/${COM_FILENAME_WITHOUT_LOC}.fastq" > "${COM_FILE_PATH}/${REF_FILENAME_WITHOUT_LOC}_${COM_FILENAME_WITHOUT_LOC}.sam"
		echo "Sam File Created with Single Read"	
	elif [ $IF_PAIRED_OR_NOT -eq 1 ]
	then
	   $2/bwa sampe -n ${11} $REF_FILENAME "$3/${REF_FILENAME_WITHOUT_LOC}_${COM_FILENAME_WITHOUT_LOC}_aln_sa1.sai" "$3/${REF_FILENAME_WITHOUT_LOC}_${COM_FILENAME_WITHOUT_LOC}_aln_sa2.sai" "${COM_FILE_PATH}/${COM_FILENAME_WITHOUT_LOC}_1.fastq" "${COM_FILE_PATH}/${COM_FILENAME_WITHOUT_LOC}_2.fastq" > "${COM_FILE_PATH}/${REF_FILENAME_WITHOUT_LOC}_${COM_FILENAME_WITHOUT_LOC}.sam"
		echo "Sam File Created with Paired End Read"
	fi

	# Delete the temporary files once the output sam file is created.
	# -f is for force removal of the file.
	if [ $IF_PAIRED_OR_NOT -eq 0 ]
	then
	   rm -f "$3/${REF_FILENAME_WITHOUT_LOC}_${COM_FILENAME_WITHOUT_LOC}_aln_sa.sai"
	elif [ $IF_PAIRED_OR_NOT -eq 1 ]
	then
	   rm -f "$3/${REF_FILENAME_WITHOUT_LOC}_${COM_FILENAME_WITHOUT_LOC}_aln_sa1.sai"
	   rm -f "$3/${REF_FILENAME_WITHOUT_LOC}_${COM_FILENAME_WITHOUT_LOC}_aln_sa2.sai"
	fi
elif [ ${10} -eq 1 ]
then
	$2/bwa bwasw $REF_FILENAME "${COM_FILE_PATH}/${COM_FILENAME_WITHOUT_LOC}.fastq" > "${COM_FILE_PATH}/${REF_FILENAME_WITHOUT_LOC}_${COM_FILENAME_WITHOUT_LOC}.sam"
	echo "Sam file created with BWASW option"
fi

# END
