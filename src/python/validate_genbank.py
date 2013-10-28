#! /usr/bin/python

"""
validate_genbank.py - Validate a genbank file
By:  Shaun Adkins (sadkins@som.umaryland.edu)

python validate_genbank.py -g /path/to/gbk.list -o /path/to/out/dir

Requires Biopython-1.62 to run

This script can be divided into 2 parts essentially:
1) Prevalidation - Things that have to be corrected so the Biopython parser will not die
2) Validation - Things that have to be corrected after parsing a Genbank file into a Biopython SeqRecord object

One can view the changes made in genbank_changelog.txt located in the specified output path.

--genbank_file, -g => Path to a Genbank-formatted file.  Genbank file must correspond to DNA sequences
--output_path, -o => Directory path to write output.  
--log_file, -l => Name of the changelog file (only the name).
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
    line_count = 0;	# Keeps track of char count at end of LOCUS line
    id_flag = 0;	# Flag to keep track of if ID was added to front of Accession ID
    base_gbk = basename(genbank)
    if base_gbk.endswith("gb"):
        base_gbk = base_gbk + "k"	# try to keep extensions uniform
    elif base_gbk.endswith("gbwithparts"):
        base_gbk = re.sub("withparts", "k", base_gbk)
    out_file = prepare + "/" + base_gbk
    genbank_h = open_file(genbank)
    genbank_h2 = open_file(genbank)	# Will keep track of the current ACCESSION line
    out_h = open(out_file, "w")
    
    # Probably could put this for loop in a function for a cleaner appearance
    for line in iter(genbank_h.readline,''):	# Iterate until EOF
        # Deal with things if LOCUS line is encountered
    	if line.startswith("LOCUS"):
    	    old_line = line
    	    line_count = genbank_h.tell()
            if re.search("dna", line):
                #log_h.write("Found 'dna' in LOCUS line and replacing with 'DNA'.\n")
                line = re.sub("dna", "DNA", line)
            elif re.search("rna", line):
                #log_h.write("Found 'rna' in LOCUS line and replacing with 'RNA'.\n")
                line = re.sub("rna", "RNA", line)	# not working on RNA but still need capitalized for parsing later
                
            if re.search("\.pseudomolecule", line):
                #log_h.write("Removing 'pseudomolecule' from locus name in LOCUS line as locus must be less than 16 characters. \n")
                line = re.sub("\.pseudomolecule", "", line)
                
            m = re.match("LOCUS\s+(\S+)\s+", line)
            id = m.group(1)
            id_tmp = id
            id = begins_with_digit(id, log_h)
            line = re.sub(id_tmp, id, line)
            
            # When we encounter a LOCUS, we need to get accession ID and modify name if starts with digit
            genbank_h2.seek(line_count)	# Start at current LOCUS line to make sure we get right corresponding ACCESSION line
            line2 = genbank_h2.readline()
            while not line2.startswith("ACCESSION"):
                line2 = genbank_h2.readline()
            m1 = re.match("ACCESSION\s+(\S+)", line2)
            if not m1:
                sys.stderr.write("ACCESSION ID is not present in Genbank file for locus line (" + line + ")\n")
                sys.exit(1)
            accession = m1.group(1)
            acc_tmp = accession
            accession = begins_with_digit(accession, log_h)
            # Did accession get "ID" written in front?
            if accession != acc_tmp:
                id_flag = 1	# Remember to change this in the line in the genbank file later
                
            # Biopython fails if locus name is longer than 16 characters    
            if len(id) > 16:	
                #log_h.write("Locus name " + id + " is longer than 16 characters... attempting to substitute with accession ID. \n")
                if len(accession) > 16:	# pointless to substitute if accession causes Biopython to fail too
                    sys.stderr.write("Cannot use Accession ID for substitution as the ID is longer than 16 chracters.  Please consult NCBI Genbank formatting standards. Locus line (" + line + ")\n")
                    sys.exit(1)
                if accession.lower() == "unknown":
                    sys.stderr.write("Cannot substitute locus name with accession ID.  Please verify your Genbank file to make sure it meets NCBI standards. Locus line (" + line + ")\n")
                    sys.exit(1)
                else:
                    #log_h.write("Replacing locus name " + id + " with accession ID " + accession + ". \n")
                    line = re.sub(id, accession, line)
            # end if
            if line != old_line:
                log_h.write("OLD LOCUS LINE: " + old_line.rstrip() + "\n")
                log_h.write("NEW LOCUS LINE: " + line.rstrip() + "\n")
    	# end if
    	
    	# Deal with things if ACCESSION line is encountered
        if line.startswith("ACCESSION") and id_flag:	# Change Accession ID in Genbank file
            log_h.write("OLD ACCESSION: " + line.rstrip() + "\n")
            line = re.sub("ACCESSION   ", "ACCESSION   ID", line)
            id_flag = 0
            log_h.write("NEW ACCESSION: " + line.rstrip() + "\n")
        # end if
        out_h.write(line)
    # end for-loop    
    
    #add other rules that would fail during Biopython parsing if needed
    genbank_h.close()
    genbank_h2.close()
    out_h.close()
    return out_file

# If string begins with a digit, add "ID" to front.  Return string regardless of change
def begins_with_digit(string, log_h):
    m = re.match("^(\d+\S+)", string)
    if m:
        string = "ID" + m.group(1)
        #log_h.write("Because it starts with a digit, " + m.group(1) + " has been replaced with " + string + ".\n")
    return string

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
        if len(gb_record.features) > 0:
            remove_genes_from_circular_starting_at_end(gb_record, log_h)
            fix_db_xref(gb_record, log_h)
        else:
            log_h.write("No annotation features present!!!  Skipping feature checks!!!")
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
        sys.stdere.write("RNA or Protein alphabet detected for Genbank file.  Must supply only nucleotide Genbank files.\n")
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
    	log_h.write("Changing DEFINITION...\nOLD VALUE: " + record.description + "\n")
        record.description = dash.sub("_", record.description)	#...substitute for an underscore
        record.description = colon.sub("_", record.description)
        record.description = comma.sub("_", record.description)
        record.description = piper.sub("_", record.description)
        log_h.write("NEW VALUE: " + record.description + "\n")
    #SOURCE        
    m2 = dash.search(record.annotations['source'])
    n2 = colon.search(record.annotations['source'])
    o2 = comma.search(record.annotations['source'])
    p2 = piper.search(record.annotations['source'])
    if m2 or n2 or o2 or p2:
    	log_h.write("Changing SOURCE...\nOLD VALUE: " + record.annotations['source'] + "\n")    
        record.annotations['source'] = dash.sub("_", record.annotations['source'])
        record.annotations['source'] = colon.sub("_", record.annotations['source'])
        record.annotations['source'] = comma.sub("_", record.annotations['source'])
        record.annotations['source'] = piper.sub("_", record.annotations['source'])        
        log_h.write("NEW VALUE: " + record.annotations['source'] + "\n")
    #ORGANISM
    m3 = dash.search(record.annotations['organism'])
    n3 = colon.search(record.annotations['organism'])
    o3 = comma.search(record.annotations['organism'])
    p3 = piper.search(record.annotations['organism'])
    if m3 or n3 or o3 or p3:
    	log_h.write("Changing ORGANISM...\nOLD VALUE: " + record.annotations['organism'] + "\n")    
        record.annotations['organism'] = dash.sub("_", record.annotations['organism'])
        record.annotations['organism'] = colon.sub("_", record.annotations['organism'])
        record.annotations['organism'] = comma.sub("_", record.annotations['organism'])
        record.annotations['organism'] = piper.sub("_", record.annotations['organism'])    
        log_h.write("NEW VALUE: " + record.annotations['organism'] + "\n")
    #FEATURES.source.organism and FEATURES.source.strain
    for feature in record.features:
        if feature.type == 'source':
            assert len(feature.qualifiers['organism']) == 1, "This record has more than one organism listed in the FEATURES.source entry"
            m4 = dash.search(feature.qualifiers['organism'][0])	# Qualifiers return as lists, but organism should only have 1 element
            n4 = colon.search(feature.qualifiers['organism'][0])
            o4 = comma.search(feature.qualifiers['organism'][0])
            p4 = piper.search(feature.qualifiers['organism'][0])
            if m4 or n4 or o4 or p4:
                log_h.write("Changing FEATURES.source.organism...\nOLD VALUE: " + feature.qualifiers['organism'][0] + "\n")                 
                feature.qualifiers['organism'][0] = dash.sub("_", feature.qualifiers['organism'][0])
                feature.qualifiers['organism'][0] = colon.sub("_", feature.qualifiers['organism'][0])
                feature.qualifiers['organism'][0] = comma.sub("_", feature.qualifiers['organism'][0])
                feature.qualifiers['organism'][0] = piper.sub("_", feature.qualifiers['organism'][0])
                log_h.write("NEW VALUE: " + feature.qualifiers['organism'][0] + "\n")                 
            
            m5 = dash.search(feature.qualifiers['strain'][0])	# Qualifiers return as lists, but strain should only have 1 element
            n5 = colon.search(feature.qualifiers['strain'][0])
            o5 = comma.search(feature.qualifiers['strain'][0])
            p5 = piper.search(feature.qualifiers['strain'][0])
            if m5 or n5 or o5 or p5:
                log_h.write("Changing FEATURES.source.strain...\nOLD VALUE: " + feature.qualifiers['strain'][0] + "\n")                 
                feature.qualifiers['strain'][0] = dash.sub("_", feature.qualifiers['strain'][0])
                feature.qualifiers['strain'][0] = colon.sub("_", feature.qualifiers['strain'][0])
                feature.qualifiers['strain'][0] = comma.sub("_", feature.qualifiers['strain'][0])
                feature.qualifiers['strain'][0] = piper.sub("_", feature.qualifiers['strain'][0])     
                log_h.write("NEW VALUE: " + feature.qualifiers['strain'][0] + "\n")                                                
            break
    return    

# Non- "AGCT" characters should be replaced with "N"
def replace_invalid_sequence_chars(record, log_h):
    seq = str(record.seq)
    assert len(seq) > 0, "No sequence present in Genbank file"
    m = re.search("[^AGCT]", seq.upper())	# Keep sequences uniform by making upper-case
    if m:
        log_h.write("Sequence has non-AGCT characters present... replacing those characters with 'N'.\n")
        seq = re.sub("[^AGCT]", "N", seq.upper())	# Possibly revise later to provide statistics of positional changes
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
            #print feature.location
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
    usage = "usage: %prog -g /path/to/file.gbk -o /path/to/out/dir"
    description = "Validate a Genbank file. Requires Biopython-1.62 to run"
    parser = OptionParser(usage=usage, description=description)
    parser.add_option("-g", "--genbank_file", help="path to a Genbank-formatted file");
    parser.add_option("-o", "--output_path", help="directory path to write output")
    parser.add_option("-l", "--log_file", help="name of the changelog file (only the name)");
    (options, args) = parser.parse_args()

    if not options.genbank_file:
        parser.error( "Genbank file path (-g) must be provided")
    if not options.output_path:
        parser.error("Output directory path (-o) must be provided")

    gbk = options.genbank_file
    
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
    if options.log_file:
        log = outdir + "/" + options.log_file
    else:
        log = outdir + "/gbk_changelog.txt"
    log_h = open(log, "w")

    gbk = gbk.rstrip()
    if not gbk.endswith("gbk") and not gbk.endswith("gb") and not gbk.endswith("gbwithparts"):	#Change into a regex later
        sys.stderr.write("File " + gbk + " does not have a proper Genbank file extension (.gbk or .gb)... skipping\n")
        sys.exit(1)
    log_h.write("Now preparing " + gbk + " for validation\n")	
    new_gbk = prevalidation(gbk, pre, log_h)
    log_h.write("Now validating " + gbk + " ...\n")
    validate_genbank(new_gbk, val, log_h)
    log_h.write("\n")	#Really wish I could send to stdout and stderr w/o writing 2 statements. 
	
    sys.stdout.write("Finished validating " + gbk + "\n")
    log_h.close()

if __name__ == '__main__':
    main()
    sys.exit(0)

