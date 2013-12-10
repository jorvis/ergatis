#!/usr/bin/env python

"""
Creates a workflow iterator for the get_bam_coverage_dx component.
"""

import argparse
import glob
import itertools
import os


def parse_CLI_arguments():
    """
    Parses any command-line arguments passed to this script and returns an
    ArgParse object.
    """
    parser = argparse.ArgumentParser('Creates a workflow iterator for the '
                                     'get_bam_coverage_dx component')
    parser.add_argument('-i', '--sam-file-list', dest='sam_file_list',
                        required=True, help='List of SAM files to run '
                        'through the get_bam_coverage_dx component.')
    parser.add_argument('-r', '--reference-db-path', dest='ref_db_path',
                        required=True, help='Folder containing BWA index '
                        'reference database(s)')
    parser.add_argument('-o', '--output-file', dest='out_file', required=True,
                        help='Desired output iterator file')

    return parser.parse_args()


def parse_input_file_list(file_list):
    """
    Parses a file list and returns a list of all files
    """
    with open(file_list) as x: files = x.readlines()
    return [x.strip() for x in files]


def get_reference_db_files(db_path):
    """
    Retrieves all reference database files from the provided path. Files should
    be in fasta format and already indexed by BWA.
    """
    ref_files = []

    for ref_file in glob.glob(os.path.join(db_path, "*.fasta")):
        ref_files.append(ref_file)

    return ref_files


def map_ref_to_sam_files(sam_files, reference_files):
    """
    Attempts to naively map reference files to SAM files based on the 
    basename for each file.
    """
    sam_to_ref_list = []
    ref_base_dict = dict([(os.path.splitext(os.path.basename(x))[0], x)
                      for x in reference_files])
    
    ## We are assuming that our input files, when we are passed more 
    ## than one reference, are in the format <FILE_NAME>.<REF BASE>.sam
    ## 
    ## So we want to map <REF BASE> to the basename of the reference files
    for sam_file in sam_files:
        sam_ext = os.path.splitext(sam_file)[0]
        sam_elts = sam_ext.split('.')
        sam_elts.remove('sorted_by_position')
        ref_file = ref_base_dict[sam_elts[-1]]

        sam_to_ref_list.append((sam_file, ref_file))
    
    return sam_to_ref_list
    

def main(args):
    sam_files = parse_input_file_list(args.sam_file_list)
    reference_files = get_reference_db_files(args.ref_db_path)

    iterator_fh = open(args.out_file, 'w')
    iterator_fh.write("\t".join(["$;I_FILE_BASE$;",
                                 "$;I_FILE_NAME$;",
                                 "$;I_FILE_PATH$;",
                                 "$;I_DATABASE_BASE$;",
                                 "$;I_DATABASE_FILE$;"]) + "\n")
    
    if len(reference_files) > 1:
        sam_ref_map = map_ref_to_sam_files(sam_files, reference_files)
    else:
        sam_ref_map = itertools.product(sam_files, db_files)

    for (sam_file, reference) in sam_ref_map:
        file_name = os.path.basename(sam_file)
        file_base = os.path.splitext(os.path.basename(sam_file))[0]
        ref_base = os.path.splitext(os.path.basename(reference))[0]

        iterator_fh.write("\t".join([file_base.replace('.sorted_by_position', ''), 
                                     file_name.replace('.sorted_by_position', ''),
                                     sam_file, 
                                     ref_base, reference]) + "\n")

    iterator_fh.close()


if __name__ == "__main__":
    main(parse_CLI_arguments())
