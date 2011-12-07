# This is a shell script to zip a few files.
# Using gzip to zip.
# Author: Ashwinkumar Ganesan.

# This is to find the path only.


find $1 -name '*M_*.sam' -exec gzip {} \;

exit
FILENAME_LEN=`expr length $1`
FILENAME_WITHOUT_LOC=`echo $1 | awk -F"/" '{print $NF}'`
FILENAME_WITHOUT_LOC_LEN=`expr length $FILENAME_WITHOUT_LOC`
END_LOC=`expr $FILENAME_LEN - $FILENAME_WITHOUT_LOC_LEN`
COM_FILE_PATH=`echo $1 | cut -c1-${END_LOC}`

echo "Location: ${COM_FILE_PATH}"

# This is to exit the script in case it is a single read file.
IF_PAIRED_OR_NOT=`ls -l ${COM_FILE_PATH}/*_1.fastq | wc -l`
if [ $IF_PAIRED_OR_NOT -eq 0 ]
then
   exit
fi


# Zipping the following files.
#gzip "${COM_FILE_PATH}/final_A2D.fna_${FILENAME_WITHOUT_LOC}.sam"
#gzip "${COM_FILE_PATH}/final_E2P.fna_${FILENAME_WITHOUT_LOC}.sam"
#gzip "${COM_FILE_PATH}/final_R2Z.fna_${FILENAME_WITHOUT_LOC}.sam"
#gzip "${COM_FILE_PATH}/hg19.fa_${FILENAME_WITHOUT_LOC}.sam"

#gzip "${COM_FILE_PATH}/final_A2D.fna_${FILENAME_WITHOUT_LOC}.sam_UM_UM.sam"
#gzip "${COM_FILE_PATH}/final_E2P.fna_${FILENAME_WITHOUT_LOC}.sam_UM_UM.sam"
#gzip "${COM_FILE_PATH}/final_R2Z.fna_${FILENAME_WITHOUT_LOC}.sam_UM_UM.sam"
#gzip "${COM_FILE_PATH}/hg19.fa_${FILENAME_WITHOUT_LOC}.sam_M_M.sam"
gzip -f ${COM_FILE_PATH}/*.sam
