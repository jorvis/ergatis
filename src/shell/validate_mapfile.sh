#!/bin/bash

source /opt/opt-packages/qiime-1.8/activate.sh

#### validate_mapfile.sh - Validates mapping file
#### Author:  Shaun Adkins - sadkins@som.umaryland.edu
#### Description:  This script both extracts barcodes from MiSeq input and assembles MiSeq reads in preparation for QIIME analyses
#### Input:     -o -> Path to output directory
####		-m -> A mapping text file with barcode information
#### Output:	A corrected mapping text file redirected to STDOUT

while getopts 'm:o:' opt
do
    case $opt in
        m)	MAPPING_FILE=$OPTARG
            ;;
        o)      OUTPUT_DIR=$OPTARG
            ;;
        ?)	printf ' usage: %s -o /path/to/output/dir/ -m /path/to/mapping.txt \n' $(basename $0) >&2
            exit 2
            ;;
    esac
done
shift $((OPTIND - 1))	# Shift out all processed options to handle any existing non-option args

# Verify FastQ list input is exactly 2 files
echo "Verifying contents of FASTQ and mapping files..."

#### Verify 3 files do exist
if [ ! -f $MAPPING_FILE ]; then
	echo "Required files not exist!! Abort" >&2
	exit; 
fi

echo "Copying mapping file to a general name..."
cp $MAPPING_FILE $OUTPUT_DIR/mapping.txt

####make sure the Mapping file is correct, warning is ok, just no errors
# Note:  validate_mapping_file.py is part of QIIME and QIIME's bin dir should be set in $PATH
echo "Validate provided mapping file using QIIME..."
validate_mapping_file.py -m $OUTPUT_DIR/mapping.txt -o $OUTPUT_DIR	## for QIIME 1.8
if (( $? )) ; then echo "Failed at validate_mapping_file.py" >&2 ; exit 1; fi

echo "Done!"
exit 0
