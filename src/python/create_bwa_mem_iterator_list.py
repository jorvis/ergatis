#!/usr/bin/env python

"""
Creates a workflow iterator for the BWA mem component.
"""

import argparse
import itertools
import os


def parse_CLI_arguments():
    """
    Parses any command-line arguments passed to this script and returns an
    ArgParse object.
    """
    parser = argparse.ArgumentParser('Creates a workflow iterator for the ' 
                                     'BWA mem component.')
    parser.add_argument('-i', '--sequence-file-list', dest='seq_file_list',
                        required=True, help='List of paired-end files to run '
                        'through BWA mem')
    parser.add_argument('-r', '--reference-file-list', dest='ref_file_list',
                        required=True, help='List of reference files to map '
                        'reads against')
    parser.add_argument('-o', '--output-file', dest='out_file', required=True,
                        help='Desired output iterator file')
                        

    return parser.parse_args()


def parse_input_files(sequence_files):
    """
    Parses list of input sequence files and returns a list containing mate 
    pairs or single-ended reads
    """
    files_list = []
    
    with open(sequence_files) as x: files = [y.strip() for y in x.readlines()]
    for line in sorted(files):
        if "_2.fastq" in line:
            if line.replace("_2.fastq", "_1.fastq") in files:
                continue
            else:
                raise Exception("ERROR: Paired-end file %s has no mate pair" % line)
        if "_1.fastq" in line:
            file_base = os.path.splitext(line)[0]
            mate_pair = line.replace("_1.fastq", "_2.fastq")

            ## Check if the mate-pair is in our raw files list
            if mate_pair in files:
                files_list.append((line, mate_pair))
            else:
                ## It looks like we have what seems to be a mate pair file 
                ## but no corresponding pair
                raise Exception("ERROR: Paired-end file %s has no mate pair" % line)

        else:
            ## Assuming we are dealing with single-ended files now
            files_list.append((line.strip(),))

    return files_list


def parse_ref_db_list(ref_file_list):
    """
    Parses a file list of reference database files (in fasta format). 
    """
    with open(ref_file_list) as x: files = x.readlines()
    return [x.strip() for x in files if x.endswith('.fasta\n')]


def main(args):
    sequence_files = parse_input_files(args.seq_file_list)
    db_files = parse_ref_db_list(args.ref_file_list)

    iterator_fh = open(args.out_file, 'w') 
    iterator_fh.write("\t".join(["$;I_FILE_BASE$;",
                                 "$;I_FILE_NAME$;",
                                 "$;I_FILE_PATH$;",
                                 "$;I_MATE_PAIR$;",
                                 "$;I_REF_BASE$;",
                                 "$;I_REF_FILE$;"]) + "\n")

    for (files, reference) in itertools.product(sequence_files, db_files):
        ## We want to check if we are dealing with mate pairs here
        if len(files) > 1:
            target_file = files[0]
            mate_pair = files[1]
        else:
            target_file = files[0]
            mate_pair = ""

        file_name = os.path.basename(target_file)
        file_base = os.path.splitext(os.path.basename(target_file))[0]
        ref_base = os.path.splitext(os.path.basename(reference))[0]
        iterator_fh.write("\t".join([file_base, file_name, target_file, 
                                     mate_pair, ref_base, reference]) + "\n")

    iterator_fh.close()

if __name__ == "__main__":
    main(parse_CLI_arguments())
