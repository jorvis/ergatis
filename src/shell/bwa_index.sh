# This program indexes the reference file using BWA.

#!/bin/sh
# AWK has been used to parse the file.
# Author: Ashwinkumar Ganesan.
# Date Created: 11-12-2010

# This is for testing purposes.
# set -vx

# Following is the list of parameters to the shell script.
# 1  - Is the File path provided by the user or a given pipeline.
# 2  - Is the location of the bwa executable file.
# 3  - Is the algorithm to be used by BWA.

# Reference File Location.
# This step with the parse the filename from the file location.
REF_NAME=`echo $1 | awk -F"/" '{print $NF}'`

# This is the whole location.
REF_FILENAME=$1
#echo $REF_FILENAME

# The next step is to create an index for the file.
# The files created will be stored in the output directory location.
$2/bwa index -p $REF_FILENAME -a $3 $REF_FILENAME

# END
