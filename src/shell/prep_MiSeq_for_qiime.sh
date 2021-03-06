#!/bin/bash

#### prep_MiSeq_for_qiime.sh - preps MiSeq data to be used in QIIME analyses
#### Author:  Shaun Adkins - sadkins@som.umaryland.edu
#### Description:  This script both extracts barcodes from MiSeq input and assembles MiSeq reads in preparation for QIIME analyses
#### Input: -f -> Two paired-end FASTQ files
####		-o -> Output directory to store script outputs
####		-b -> Directory path for the $BIN_DIR/fq_mergelines.pl, $BIN_DIR/fq_splitlines.pl, and $BIN_DIR/fq_getPairAndOrphan1.8.py scripts
#### Note:  This is compatible with 1-step PCR MiSeq data, but not 2-step PCR

#### This pipeline is to process dual barcoding reads separately
#### It depends on the following package/module
#### 1. FASTX-Toolkit from http://hannonlab.cshl.edu/fastx_toolkit/
#### 2. seqtk from https://github.com/lh3/seqtk
#### 3. prinseq-lite from http://prinseq.sourceforge.net/

#### Originally was a script written by Bing Ma, but I modified to work in the CloVR 16S pipeline in Ergatis

while getopts 'b:f:o:' opt
do
    case $opt in
        b)	BIN_DIR=$OPTARG
            ;;
        f)	FASTQ_LIST=$OPTARG
            ;;
        o)	OUTPUT_DIR=$OPTARG
            ;;
        ?)	printf ' usage: %s -f /path/to/fastq.list -o /path/to/output/dir/ -b /dir/to/misc/scripts/ \n' $(basename $0) >&2
            exit 2
            ;;
    esac
done
shift $((OPTIND - 1))	# Shift out all processed options to handle any existing non-option args

# Verify FastQ list input is exactly 2 files
echo "Verifying contents of FASTQ list ..."
line_count=$(wc -l $FASTQ_LIST | cut -f1 -d' ')
if [ $line_count -ne 2 ]; then
	echo "FASTQ list must have exactly 2 sequences in it.  This file had $line_count" >&2
	exit 1
fi

R1=`head -n 1 $FASTQ_LIST`
R2=`tail -n 1 $FASTQ_LIST`


#### Verify 3 files do exist
if [ ! -f $R1 ] || [ ! -f $R2 ] ; then
	echo "Required files not exist!! Abort" >&2
	exit; 
fi

length_with_barcode=250 ## for 250PE, this variable is not used in the script, just for reference
### cutoffLength is about 0.75 of the length_with_barcode
cutoffLength=187 ### ie, 250 x 0.75

prep_dir=$OUTPUT_DIR/prep

if [ ! -d $prep_dir ]; then
	mkdir -m 777 $prep_dir;
fi

#trim barcode, form a new barcode file 
echo "Extracting barcodes from reads from FASTQ input to create barcode FASTQ files..."
fastx_trimmer -i $R1 -f 1 -l 12 -Q 33 -o $prep_dir/R1_barcode.fq
fastx_trimmer -i $R2 -f 1 -l 12 -Q 33 -o $prep_dir/R2_barcode.fq
if (( $? )) ; then echo "Failed at fastx-trimmer"; exit 1 ; fi

echo "Creating tab-delimited files from the barcode FASTQ files..."
cat $prep_dir/R1_barcode.fq | $BIN_DIR/fq_mergelines.pl > $prep_dir/R1_barcode_temp
cat $prep_dir/R2_barcode.fq | $BIN_DIR/fq_mergelines.pl > $prep_dir/R2_barcode_temp

echo "Merging barcodes of both paired-end files..."
paste $prep_dir/R1_barcode_temp $prep_dir/R2_barcode_temp | awk -F"\t" '{print $5"\t"$2$6"\t"$3"\t"$4$8}' | $BIN_DIR/fq_splitlines.pl > $prep_dir/R1R2_barcode.fastq
cat $prep_dir/R1R2_barcode.fastq | $BIN_DIR/fq_mergelines.pl > $prep_dir/R1R2_barcode_temp

#trim sequences (12 bp of barcode) off the sequence
echo "Trimming off barcodes from FASTQ input..." 
seqtk trimfq -b 12 $R1 > $prep_dir/R1_trimmed_seq.fastq
seqtk trimfq -b 12 $R2 > $prep_dir/R2_trimmed_seq.fastq
if (( $? )) ; then echo "Failed at seqtk barcode trimming"; exit 1 ; fi

#trim off low quality bases 
echo "Trimming off low quality bases (error threshold = 0.05)..."
seqtk trimfq $prep_dir/R1_trimmed_seq.fastq | $BIN_DIR/fq_mergelines.pl > $prep_dir/R1_seqtk_trimmed_seq_temp  
seqtk trimfq $prep_dir/R2_trimmed_seq.fastq | $BIN_DIR/fq_mergelines.pl > $prep_dir/R2_seqtk_trimmed_seq_temp 

#create new fastq file to concatenate barcode (24bp) with sequences
echo "Creating new FASTQ inputs with concatenated barcodes and input sequences..."
paste $prep_dir/R1R2_barcode_temp $prep_dir/R2_seqtk_trimmed_seq_temp | awk -F"\t" '{print $1"\t"$2$6"\t"$7"\t"$4$8}' | $BIN_DIR/fq_splitlines.pl | sed 's/2:N:0:/3:N:0:/' > $prep_dir/R3N.fq
paste $prep_dir/R1R2_barcode_temp $prep_dir/R1_seqtk_trimmed_seq_temp | awk -F"\t" '{print $1"\t"$2$6"\t"$7"\t"$4$8}' | $BIN_DIR/fq_splitlines.pl | sed 's/2:N:0:/1:N:0:/' > $prep_dir/R1N.fq

####based on the length, define the cutoff of the concat length, default is half of the original 250bp
echo "Running prinseq-lite to filter out bad reads..."
prinseq-lite.pl -fastq $prep_dir/R3N.fq -out_format 3 -out_good $prep_dir/R3N_${cutoffLength} -no_qual_header -line_width 0 -min_len $cutoffLength
prinseq-lite.pl -fastq $prep_dir/R1N.fq -out_format 3 -out_good $prep_dir/R1N_${cutoffLength} -no_qual_header -line_width 0 -min_len $cutoffLength
if (( $? )) ; then echo "Failed at prinseq-lite.pl"; exit 1 ; fi

###in this case, the total and total_filter files have no difference, only the barcode file will change. The purpose is to make the trimmed seq file and barcode file match. 
# Note - Requires BioPython v 1.51 or later
echo "Creating paired-end and barcode FASTQ files"
$BIN_DIR/fq_getPairAndOrphan.py $prep_dir/R3N_${cutoffLength}.fastq $prep_dir/R1N_${cutoffLength}.fastq $OUTPUT_DIR/R3N_PE.fq $OUTPUT_DIR/R1N_PE.fq $prep_dir/orphan_temp.fq 
$BIN_DIR/fq_getPairAndOrphan.py $OUTPUT_DIR/R3N_PE.fq $prep_dir/R1R2_barcode.fastq $prep_dir/PE_temp.fq $OUTPUT_DIR/barcode.fq $prep_dir/orphan_temp.fq 
if (( $? )) ; then echo "Failed at fq_getPairAndOrphan.py"; exit 1 ; fi

# The files we want to pass into the QIIME (1.8) split_libraries_fastq.py script are
# 1) R3N_PE.fq
# 2) R1N_PE.fq
# 3) barcode.fq

echo "Removing temporary files..."
rm $prep_dir/*_bad_*
rm $prep_dir/*_temp*
#rm -rf $prep_dir
#rm R1_barcode.fq R2_barcode.fq R1R2_barcode.fastq R1_trimmed_seq.fastq R2_trimmed_seq.fastq R3N.fq R1N.fq R3N_${cutoffLength}.fastq R1N_${cutoffLength}.fastq

echo "Done!"
exit 0
