#! /usr/bin/python

"""
validate_genbank.py - Validate a list of genbank files (must be nucleotide files)
By:  Shaun Adkins (sadkins@som.umaryland.edu)

python validate_genbank.py -g /path/to/gbk.list -o /path/to/out/dir

--genbank_list, -g => A line-delimited list of Genbank file paths.  Genbank files must correspond to nucleotide sequences
--output_path, -o => Directory path to write output
"""

import sys
import os
import re
from optparse import OptionParser	#may switch to argparse in the future
from os.path import basename
from Bio import SeqIO
from Bio import SeqFeature
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

#############
#      FUNCTIONS      #
#############

# Some files need to have changes made so they won't fail during Biopython parsing. 
# This is where we fix that
def prevalidation(genbank, prepare, log_h):
    base_gbk = basename(genbank)
    if base_gbk.endswith("gb"):
        base_gbk = base_gbk + "k"	# try to keep extensions uniform
    elif base_gbk.endswith("gbwithparts"):
        base_gbk = re.sub("withparts", "k", base_gbk)
    out_file = prepare + "/" + base_gbk
    genbank_h = open_file(genbank)
    out_h = open(out_file, "w")
    for line in genbank_h:
    	if line.startswith("LOCUS"):
            if re.search("dna", line):
                log_h.write("Found 'dna' in LOCUS line and replacing with 'DNA'.\n")
                line = re.sub("dna", "DNA", line)
            elif re.search("rna", line):
                log_h.write("Found 'rna' in LOCUS line and replacing with 'RNA'.\n")
                line = re.sub("rna", "RNA", line)	# not working on RNA but still need capitalized for parsing later
            if re.search("\.pseudomolecule", line):
                log_h.write("Removing 'pseudomolecule from locus name in LOCUS line as locus must be less than 16 characters. \n")
                line = re.sub("\.pseudomolecule", "", line)
        out_h.write(line)
    #add other rules that would fail during Biopython parsing if needed
    genbank_h.close()
    out_h.close()
    return out_file

# Using Biopython to validate genbank files
def validate_genbank(genbank, valid, log_h):
    base_gbk = basename(genbank)
    print base_gbk
    out_f = valid + "/" + base_gbk
    genbank_h = open_file(genbank)
    record_list = parse_file(genbank_h)
    for gb_record in record_list:
        is_sequence_nucleotide(gb_record, log_h)
        is_accession_present(gb_record, log_h)
        replace_invalid_header_chars(gb_record, log_h)
        replace_invalid_sequence_chars(gb_record, log_h)
        remove_genes_from_circular_starting_at_end(gb_record, log_h)
        fix_db_xref(gb_record, log_h)
        log_h.write("\n")	# newline after contigs
    #add other rules as we expand this script
    write_output(record_list, out_f)

# Opens genbank file
def open_file(gb):
    try:
        gb_h = open(gb, "r")
    except IOError:
        sys.stderr.write("Cannot open file " + gb + " for reading!\n")
        sys.exit(1)
    else:
        return gb_h

# Parses the given Genbank file and creates a SeqRecord object of it
def parse_file(gb_h):
    record_list = []
    for record in SeqIO.parse(gb_h, "genbank"):	# parse genbank into a SeqRecord object
        record_list.append(record)
        #print record
    gb_h.close()
    return record_list

# Checks to make sure the sequence alphabet is DNA
def is_sequence_nucleotide(record, log_h):
    #print str(record.seq.alphabet)
    m = re.search("RNA", str(record.seq.alphabet))
    n = re.search("Protein", str(record.seq.alphabet))
    if m or n:
        log_h.write("RNA or Protein alphabet detected for Genbank file.  Must supply only nucleotide Genbank files.\n")
        sys.exit(1)

# Accession IDs should be present in every Genbank file and meet proper format
def is_accession_present(record, log_h):
    try:
        log_h.write("Accession ID: " + record.id + "\n")
    except NameError:
        log_h.write("Accession ID not found!!!\n")
    """else:
        p1 = re.compile("[a-zA-Z]{1}_?\d{5}")	# two separate types of nucleotide Accession IDs
        p2 = re.compile("[a-zA-Z]{2}_?\d{6}")
        m1 = p1.match(record.id)
        m2 = p2.match(record.id)
        if m1 or m2:
            log_h.write("Valid nucleotide accession ID\n")
        else:
            log_h.write("Accession ID: " + record.id + " is not valid.  A nucleotide-based accession ID from Genbank must have 2 letters and 6 digits (LL######) or 1 letter and 5 digits (L#####).  RefSeq accession IDs have an underscore in the 3rd position\n")
    """
    return

# Organism name and Features, Source, and Organism attributes need "-" or ":" replaced with "_"
def replace_invalid_header_chars(record, log_h):
    dash = re.compile("-")	#compiling patterns for both the dash and the colon
    colon = re.compile(":")
    comma = re.compile(",")
    piper = re.compile("\|")	#naming piper in case 'pipe' is a key word
    
    #DEFINITION    
    m1 = dash.search(record.description)	#searching for matches
    n1 = colon.search(record.description)
    o1 = comma.search(record.description)
    p1 = piper.search(record.description)
    if m1 or n1 or o1 or p1:	#if match was found for either dash or colon...
        log_h.write("A dash (-), colon (:), pipe (|), or comma (,) is present in the DEFINITION section and will be converted into an underscore (_).\n")
        record.description = dash.sub("_", record.description)	#...substitute for an underscore
        record.description = colon.sub("_", record.description)
        record.description = comma.sub("_", record.description)
        record.description = piper.sub("_", record.description)
    #SOURCE        
    m2 = dash.search(record.annotations['source'])
    n2 = colon.search(record.annotations['source'])
    o2 = comma.search(record.annotations['source'])
    p2 = piper.search(record.annotations['source'])
    if m2 or n2 or o2 or p2:
        log_h.write("A dash (-), colon (:), pipe (|), or comma (,) is present in the SOURCE section and will be converted into an underscore (_).\n")
        record.annotations['source'] = dash.sub("_", record.annotations['source'])
        record.annotations['source'] = colon.sub("_", record.annotations['source'])
        record.annotations['source'] = comma.sub("_", record.annotations['source'])
        record.annotations['source'] = piper.sub("_", record.annotations['source'])        
    #ORGANISM
    m3 = dash.search(record.annotations['organism'])
    n3 = colon.search(record.annotations['organism'])
    o3 = comma.search(record.annotations['organism'])
    p3 = piper.search(record.annotations['organism'])
    if m3 or n3 or o3 or p3:
        log_h.write("A dash (-), colon (:), pipe (|), or comma (,) is present  in the ORGANISM section and will be converted into an underscore (_).\n")
        record.annotations['organism'] = dash.sub("_", record.annotations['organism'])
        record.annotations['organism'] = colon.sub("_", record.annotations['organism'])
        record.annotations['organism'] = comma.sub("_", record.annotations['organism'])
        record.annotations['organism'] = piper.sub("_", record.annotations['organism'])    
    #FEATURES.source.organism
    for feature in record.features:
        if feature.type == 'source':
            assert len(feature.qualifiers['organism']) == 1, "This record has more than one organism listed in the FEATURES.source entry"
            m4 = dash.search(feature.qualifiers['organism'][0])	# Qualifiers return as lists, but organism should only have 1 element
            n4 = colon.search(feature.qualifiers['organism'][0])
            o4 = comma.search(feature.qualifiers['organism'][0])
            p4 = piper.search(feature.qualifiers['organism'][0])
            if m4 or n4 or o4 or p4:
                log_h.write("A dash (-), colon (:), pipe (|), or comma (,) is present in the FEATURES.source.organism section and will be converted into an underscore (_).\n")
                feature.qualifiers['organism'][0] = dash.sub("_", feature.qualifiers['organism'][0])
                feature.qualifiers['organism'][0] = colon.sub("_", feature.qualifiers['organism'][0])
                feature.qualifiers['organism'][0] = comma.sub("_", feature.qualifiers['organism'][0])
                feature.qualifiers['organism'][0] = piper.sub("_", feature.qualifiers['organism'][0])
            break
    return    

# Non- "AGCT" characters should be replaced with "N"
def replace_invalid_sequence_chars(record, log_h):
    seq = str(record.seq)
    assert len(seq) > 0, "No sequence present in Genbank file"
    m = re.search("[^AGCT]", seq.upper())	# Keep sequences uniform by making upper-case
    if m:
        log_h.write("Sequence has non-AGCT characters present... replacing those characters with 'N'.\n")
        seq = re.sub("[^AGCT]", "N", seq.upper())
        record.seq = Seq(seq.lower(), IUPAC.ambiguous_dna)	# Believe these are parsed by SeqIO as IUPACAmbiguousDNA alphabets
    #print record.seq
    #print record.seq.alphabet
    return
    
# If the genbank file has joined DNA coordinates that start at the end of a sequence
# and continue at the beginning (in circular DNA) then remove that SeqFeature
def remove_genes_from_circular_starting_at_end(record, log_h):
    f_list = []
    for feature in record.features:
        if feature.location_operator == 'join':	# Skip non-joined sequences
            flag = 0
            print feature.location
            #print len(feature.location)
            for i in range(len(feature.location)):	# Iterate through all coordinates of the list
                if i+1 < len(feature.location):	# Do not let last index run out of bounds
                    if feature.location.strand == 0:
                        if list(feature.location)[i] > list(feature.location)[i+1]:	# if prev coordinate is larger than next coord
                            #print feature
                            flag = 1
                            if 'locus_tag' in feature.qualifiers:
                                log_h.write("Gene feature with locus_tag '" + feature.qualifiers['locus_tag'][0] + 
    	                        "' has coordinates that run from the end of the circular DNA back to the beginning.  Deleting feature since this may cause issues later on.\n")
    	                    else:	# could not find locus_tag (i.e. misc features section)
    	                        log_h.write("Gene feature with coordinates " + str(feature.location) +
    	                        " runs from end of the circular DNA back to the beginning.  Deleting feature since this may cause issues later on.\n")
                            break
                    else:	# Handle complementary strands
                        if list(feature.location)[i] < list(feature.location)[i+1]:	# if prev coordinate is smaller than next coord
                            #print feature
                            flag = 1
                            if 'locus_tag' in feature.qualifiers:
                                log_h.write("Gene feature with locus_tag '" + feature.qualifiers['locus_tag'][0] + 
                                "' has coordinates that run from the end of the circular DNA back to the beginning.  Deleting feature since this may cause issues later on.\n")
                            else:
    	                        log_h.write("Gene feature with coordinates " + str(feature.location) +
    	                        " runs from end of the circular DNA back to the beginning.  Deleting feature since this may cause issues later on.\n")                                
                            break
            if flag == 0:
                f_list.append(feature)
        else:
            f_list.append(feature)
    record.features = f_list	# assign updated feature list to the record
    return
    
# If db_xref is not in the form of 'database:identifier', change it to be that.    
def fix_db_xref(record, log_h):	
    for feature in record.features: # a list of SeqFeature objects
    	if 'db_xref' not in feature.qualifiers:
    	    continue
    	db = []	# db_xref qualifier can have multiple elements so need a list
        for dbxref in feature.qualifiers['db_xref']:
            m = re.search(":", dbxref)
            if not m:
                log_h.write("DB_xref entry '" + dbxref + "' is invalid.  Adding database identifier to make it valid.\n")
                dbxref = "UNKNOWN:" + dbxref
        	dbxref = re.sub("\s+.*", '', dbxref)	# removing any information after the first space.  This may need to be modified later
        	#print dbxref
            db.append(dbxref)
        feature.qualifiers['db_xref'] = db
        #print feature.qualifiers['db_xref']
    return    
    
# Using our up-to-date Genbank record information to write a new Genbank file
def write_output(record_list, outfile):
    out_h = open(outfile, "w")
    SeqIO.write(record_list, out_h, "genbank")
    out_h.close()
    return

#######
#   MAIN   #
#######

def main():
    # Set up options parser and help usage statement
    usage = "usage: %prog -g /path/to/gbk.list -o /path/to/out/dir"
    description = "Validate a list of Genbank input files"
    parser = OptionParser(usage=usage, description=description)
    parser.add_option("-g", "--genbank_list", help="a line-delimited list of Genbank file paths");
    parser.add_option("-o", "--output_path", help="directory path to write output")
    (options, args) = parser.parse_args()

    if not options.genbank_list:
        parser.error( "Genbank list path (-g) must be provided")
    if not options.output_path:
        parser.error("Output directory path (-o) must be provided")

    f = open(options.genbank_list, "r")	#open the genbank_list file for reading
    lines = f.readlines()
    assert len(lines) > 0, "Genbank list contains no contents!"

    # Create output directory if it doesn't exist.
    outdir = options.output_path
    if not os.path.exists(outdir):
        os.mkdir(outdir, 0777)
    pre = outdir + "/prevalidate"
    val = outdir + "/validate"
    if not os.path.exists(pre):	# Directory to store modified gbk files cleaned up before file validation
        os.mkdir(pre, 0777)
    if not os.path.exists(val):
        os.mkdir(val, 0777)
    log = outdir + "/changelog.txt"
    log_h = open(log, "w")

    for gbk in lines:
        gbk = gbk.rstrip()
        if not gbk.endswith("gbk") and not gbk.endswith("gb") and not gbk.endswith("gbwithparts"):	#Change into a regex later
            log_h.write("File " + gbk + " does not have a proper Genbank file extension (.gbk or .gb)... skipping\n")
            continue
        log_h.write("Now preparing " + gbk + " for validation\n")	
        new_gbk = prevalidation(gbk, pre, log_h)
        log_h.write("Now validating " + gbk + " ...\n")
        validate_genbank(new_gbk, val, log_h)
        log_h.write("\n")	#Really wish I could send to stdout and stderr w/o writing 2 statements. 
	
    sys.stdout.write("Finished validating!\n")
    f.close()
    log_h.close()

if __name__ == '__main__':
    main()
    sys.exit(0)

