# This is shell script to delete the FASTQ files which are generated.
# Author: Ashwinkumar Ganesan.

# forcing the removal of any fastq or sam files that aren't compressed.
find $1 -name '*.fastq' -exec rm {} \;
find $1 -name '*.sam' -exec rm {} \;
find $1 -name '*.sai' -exec rm {} \;





exit
FILENAME_LEN=`expr length $1`
FILENAME_WITHOUT_LOC=`echo $1 | awk -F"/" '{print $NF}'`
FILENAME_WITHOUT_LOC_LEN=`expr length $FILENAME_WITHOUT_LOC`
END_LOC=`expr $FILENAME_LEN - $FILENAME_WITHOUT_LOC_LEN`
COM_FILE_PATH=`echo $1 | cut -c1-${END_LOC}`

# This is to exit the script in case it is a single read file.
IF_PAIRED_OR_NOT=`ls -l ${COM_FILE_PATH}/*_1.fastq | wc -l`

# Delete the fastq files at the location.
rm -f ${COM_FILE_PATH}/*.fastq

if [ $IF_PAIRED_OR_NOT -eq 0 ]
then
   touch ${COM_FILE_PATH}/single-read-file.txt
   exit
fi

echo "The Deletion is Done!!"

# End of program.
