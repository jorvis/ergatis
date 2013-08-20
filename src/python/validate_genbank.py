#! /usr/local/packages/Python-2.7/bin/python

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
from optparse import OptionParser	#if we upgrade Python beyond v2.7 switch to argparse
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
def prevalidation(genbank, prepare):
    base_gbk = basename(genbank)
    out_file = prepare + "/" + base_gbk
    genbank_h = open_file(genbank)
    out_h = open(out_file, "w")
    for line in genbank_h:
    	if line.startswith("LOCUS"):
    	    line = re.sub("dna", "DNA", line)
    	    line = re.sub("rna", "RNA", line)	# not working on RNA but still need capitalized for parsing later
    	    line = re.sub("\.pseudomolecule", "", line)
    	out_h.write(line)
    #add other rules that would fail during Biopython parsing if needed
    genbank_h.close()
    out_h.close()
    return out_file

# Using Biopython to validate genbank files
def validate_genbank(genbank, outdir):
    base_gbk = basename(genbank)
    out_f = outdir + "/" + base_gbk
    genbank_h = open_file(genbank)
    record_list = parse_file(genbank_h)
    for gb_record in record_list:
        is_sequence_nucleotide(gb_record)
        is_accession_present(gb_record)
        replace_invalid_header_chars(gb_record)
        replace_invalid_sequence_chars(gb_record)
        fix_db_xref(gb_record)
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
def is_sequence_nucleotide(record):
    #print str(record.seq.alphabet)
    m = re.search("RNA", str(record.seq.alphabet))
    n = re.search("Protein", str(record.seq.alphabet))
    if m or n:
        sys.stderr.write("RNA or Protein alphabet detected for Genbank file.  Must supply only nucleotide Genbank files.\n")
        sys.exit(1)

# Accession IDs should be present in every Genbank file and meet proper format
def is_accession_present(record):
    try:
        sys.stdout.write("Accession ID: " + record.id + "\n")
    except NameError:
        sys.stderr.write("Accession ID not found!!!\n")
    """else:
        p1 = re.compile("[a-zA-Z]{1}_?\d{5}")	# two separate types of nucleotide Accession IDs
        p2 = re.compile("[a-zA-Z]{2}_?\d{6}")
        m1 = p1.match(record.id)
        m2 = p2.match(record.id)
        if m1 or m2:
            sys.stdout.write("Valid nucleotide accession ID\n")
        else:
            sys.stderr.write("Accession ID: " + record.id + " is not valid.  A nucleotide-based accession ID from Genbank must have 2 letters and 6 digits (LL######) or 1 letter and 5 digits (L#####).  RefSeq accession IDs have an underscore in the 3rd position\n")
    """
    return

# Organism name and Features, Source, and Organism attributes need "-" or ":" replaced with "_"
def replace_invalid_header_chars(record):
    dash = re.compile("-")	#compiling patterns for both the dash and the colon
    colon = re.compile(":")
    
    #DEFINITION    
    m1 = dash.search(record.description)	#searching for matches
    n1 = colon.search(record.description)
    if m1 or n1:	#if match was found for either dash or colon...
        sys.stderr.write("A dash(-) and/or colon(:) is present in the DEFINITION section and will be converted into an underscore (_).\n")
        record.description = dash.sub("_", record.description)	#...substitute for an underscore
        record.description = colon.sub("_", record.description)
    #SOURCE        
    m2 = dash.search(record.annotations['source'])
    n2 = colon.search(record.annotations['source'])
    if m2 or n2:
        sys.stderr.write("A dash(-) and/or colon(:) is present in the SOURCE section and will be converted into an underscore (_).\n")
        record.annotations['source'] = dash.sub("_", record.annotations['source'])
        record.annotations['source'] = colon.sub("_", record.annotations['source'])
    #ORGANISM
    m3 = dash.search(record.annotations['organism'])
    n3 = colon.search(record.annotations['organism'])
    if m3 or n3:
        sys.stderr.write("A dash(-) and/or colon(:) is present in the ORGANISM section and will be converted into an underscore (_).\n")
        record.annotations['organism'] = dash.sub("_", record.annotations['organism'])
        record.annotations['organism'] = colon.sub("_", record.annotations['organism'])    
    #FEATURES.source.organism
    for feature in record.features:
        if feature.type == 'source':
            assert (len(feature.qualifiers['organism']) == 0, "This record has more than one organism listed in the FEATURES.source entry")
            m4 = dash.search(feature.qualifiers['organism'][0])	# Qualifiers return as lists, but organism should only have 1 element
            n4 = colon.search(feature.qualifiers['organism'][0])
            if m4 or n4:
                sys.stderr.write("A dash(-) and/or colon(:) is present in the FEATURES.source.organism section and will be converted into an underscore (_).\n")
                feature.qualifiers['organism'][0] = dash.sub("_", feature.qualifiers['organism'][0])
                feature.qualifiers['organism'][0] = colon.sub("_", feature.qualifiers['organism'][0])
            break
    return    

# Non- "AGCT" characters should be replaced with "N"
def replace_invalid_sequence_chars(record):
    seq = str(record.seq)
    assert (len(seq) == 0, "No sequence present in Genbank file")
    m = re.search("[^AGCT]", seq.upper())	# Keep sequences uniform by making upper-case
    if m:
        sys.stderr.write("Sequence has non-AGCT characters present... replacing those characters with 'N'.\n")
        seq = re.sub("[^AGCT]", "N", seq.upper())
        record.seq = Seq(seq.lower(), IUPAC.ambiguous_dna)	# Believe these are parsed by SeqIO as IUPACAmbiguousDNA alphabets
    #print record.seq
    #print record.seq.alphabet
    return
    
# If db_xref is not in the form of 'database:identifier', change it to be that.    
def fix_db_xref(record):	
    for feature in record.features: # a list of SeqFeature objects
    	if 'db_xref' not in feature.qualifiers:
    	    continue
    	db = []	# db_xref qualifier can have multiple elements so need a list
        for dbxref in feature.qualifiers['db_xref']:
            m = re.search(":", dbxref)
            if not m:
                sys.stderr.write("DB_xref entry " + dbxref + " is invalid.  Adding database identifier to make it valid.\n")
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
    for record in record_list:
        SeqIO.write(record, out_h, "genbank")
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
    assert (len(lines) > 0, "Genbank list contains no contents!")

    # Create output directory if it doesn't exist.
    outdir = options.output_path
    if not os.path.exists(outdir):
        os.mkdir(outdir, 0777)
    pre = outdir + "/prevalidate"
    if not os.path.exists(pre):	# Directory to store modified gbk files cleaned up before file validation
        os.mkdir(pre, 0777)

    for gbk in lines:
        gbk = gbk.rstrip()
        if not gbk.endswith("gbk") and not gbk.endswith("gb"):	#Change into a regex later
            sys.stderr.write("File " + gbk + " does not have a proper Genbank file extension (.gbk or .gb)... skipping\n")
            continue
        sys.stdout.write("Now preparing " + gbk + " for validation\n")	
        sys.stderr.write("Now preparing " + gbk + " for validation\n")
        new_gbk = prevalidation(gbk, pre)
        sys.stdout.write("Now validating " + gbk + " ...\n")
        sys.stderr.write("Now validating " + gbk + "...\n")
        validate_genbank(new_gbk, outdir)
        sys.stdout.write("\n")	#Really wish I could send to stdout and stderr w/o writing 2 statements. 
        sys.stderr.write("\n") 	#Perhaps I'll utilize the 'logging' module if I do an update.
	
    sys.stdout.write("Finished validating!\n")
    f.close()

if __name__ == '__main__':
    main()
    sys.exit(0)

