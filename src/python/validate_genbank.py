#! /usr/bin/python

"""
validate_genbank.py - Validate a genbank file
By:  Shaun Adkins (sadkins@som.umaryland.edu)

python validate_genbank.py -g /path/to/file.gbk -o /path/to/out/dir

Requires Biopython-1.62 to run

This script can be divided into 2 parts essentially:
1) Prevalidation - Things that have to be corrected so the Biopython parser will not die
2) Validation - Things that have to be corrected after parsing a Genbank file into a Biopython SeqRecord object

One can view the changes made in genbank_changelog.txt located in the specified output path.

--genbank_file, -g => Path to a Genbank-formatted file.  Genbank file must correspond to DNA sequences
--output_path, -o => Directory path to write output.
--log_file, -l => Name of the changelog file (only the name).
--debug, -d => Run in debug mode
"""

import sys
import os
import re
from optparse import OptionParser	#may switch to argparse in the future
from os.path import basename
from os.path import splitext
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
    base_gbk = re.sub("\.gb\w*",".gbk", base_gbk)	# Keeping various extensions uniform

    if os.stat(genbank).st_size == 0:	# Die if genbank file is empty
        sys.stderr.write("Genbank file " + genbank + " is of 0 size\n")
        sys.exit(1)
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
                line = re.sub("rna", "RNA", line)	# scope of program is not RNA but still need capitalized for parsing later

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
                sys.stderr.write("ACCESSION ID is not present in Genbank file for locus line (" + line + ")\nFile: " + genbank + "\n")
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
                    sys.stderr.write("Cannot use Accession ID for substitution as the ID is longer than 16 chracters.  Please consult NCBI Genbank formatting standards. Locus line (" + line + ")\nFile: " + genbank + "\n")
                    sys.exit(1)
                if accession.lower() == "unknown":
                    sys.stderr.write("Cannot substitute locus name with accession ID.  Please verify your Genbank file to make sure it meets NCBI standards. Locus line (" + line + ")\nFile: " + genbank + "\n")
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
            log_h.write("NEW ACCESSION: " + line.rstrip() + "\n")
        # end if

        if line.startswith("VERSION") and id_flag:	# Must change Version ID to match Accession
            log_h.write("OLD VERSION: " + line.rstrip() + "\n")
            line = re.sub("VERSION     ", "VERSION     ID", line)
            id_flag = 0
            log_h.write("NEW VERSION: " + line.rstrip() + "\n")

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
    organism = {}	# Initialize empty dictionary of organisms
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
            are_gene_features_present(gb_record, log_h)
#            remove_trna_and_rrna_features(gb_record, log_h)
            remove_genes_from_circular_starting_at_end(gb_record, log_h)
            fix_db_xref(gb_record, log_h)
        else:
            log_h.write("No annotation features present!!!  Skipping feature checks!!!")
        log_h.write("\n")	# newline after contigs
        new_org, trunc_org = format_organism_name(gb_record.annotations['organism'], log_h)
        # Add edited organism information to dictionary as [new_org, trunc_org]
        organism.setdefault(gb_record.annotations['organism'], []).append(new_org)
        organism.setdefault(gb_record.annotations['organism'], []).append(trunc_org)
    #add other rules as we expand this script
    write_output(record_list, out_f)
    overwrite_organism(out_f, organism)

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
    pattern = re.compile("[^A-Za-z0-9_\s]")
    m1 = pattern.search(record.description)	#searching for pattern matches
    m2 = pattern.search(record.annotations['source'])
    m3 = pattern.search(record.annotations['organism'])
    m4 = pattern.search(record.name)

    if m1:	# If match was found...substitute pattern with an underscore
        log_h.write("Changing DEFINITION...\nOLD VALUE: " + record.description + "\n")
        record.description = pattern.sub("_", record.description)
        log_h.write("NEW VALUE: " + record.description + "\n")
    if m2:
        log_h.write("Changing SOURCE...\nOLD VALUE: " + record.annotations['source'] + "\n")
        record.annotations['source'] = pattern.sub("_", record.annotations['source'])
        log_h.write("NEW VALUE: " + record.annotations['source'] + "\n")
    if m3:
        log_h.write("Changing ORGANISM...\nOLD VALUE: " + record.annotations['organism'] + "\n")
        record.annotations['organism'] = pattern.sub("_", record.annotations['organism'])
        log_h.write("NEW VALUE: " + record.annotations['organism'] + "\n")
    if m4:
        log_h.write("Changing LOCUS name...\nOLD VALUE: " + record.name + "\n")
        record.name = pattern.sub("_", record.name)
        log_h.write("NEW VALUE: " + record.name + "\n")

    # Deal with pattern changes in SeqRecord features
    for feature in record.features:
        if feature.type == 'source':
            assert len(feature.qualifiers['organism']) == 1, "This record has more than one organism listed in the FEATURES.source entry"
            if 'organism' in feature.qualifiers:
                # Qualifiers return as lists, but organism/strain should only have 1 element
                m6 = pattern.search(feature.qualifiers['organism'][0])
                if m6:
                    log_h.write("Changing FEATURES.source.organism...\nOLD VALUE: " + feature.qualifiers['organism'][0] + "\n")
                    feature.qualifiers['organism'][0] = pattern.sub("_", feature.qualifiers['organism'][0])
                    log_h.write("NEW VALUE: " + feature.qualifiers['organism'][0] + "\n")
            if 'strain' in feature.qualifiers:
                m7 = pattern.search(feature.qualifiers['strain'][0])
     	        if m7:
                    log_h.write("Changing FEATURES.source.strain...\nOLD VALUE: " + feature.qualifiers['strain'][0] + "\n")
                    feature.qualifiers['strain'][0] = pattern.sub("_", feature.qualifiers['strain'][0])
                    log_h.write("NEW VALUE: " + feature.qualifiers['strain'][0] + "\n")
        break	# Break 'for' loop after source-type features are read

    return

# Non- "AGCT" characters should be replaced with "N"
def replace_invalid_sequence_chars(record, log_h):
    seq = str(record.seq)

	# If no sequence is present, the BioPython parser will create a sequence of N's to take its place
    m1 = re.search("[AGCT]", seq.upper())
    if not m1:
    	sys.stderr.write("No ORIGIN sequence present for " + record.name + "\n")
    	sys.exit(1)

    # Replace all non-AGCT chars with N if any exist.
    p = re.compile("[^AGCT]")
    count = p.findall(seq.upper())
    if len(count) > 0:
        log_h.write("Found " + str(len(count)) + " instance(s) of non-ACGT characters present in the sequence...replaceing those characters with 'N'.\n")
        seq = p.sub("N", seq.upper())
        record.seq = Seq(seq.lower(), IUPAC.ambiguous_dna)	# Believe these are parsed by SeqIO as IUPACAmbiguousDNA alphabets

    #print record.seq
    #print record.seq.alphabet
    return

# Checks for presence of at least one gene feature annotation.  Will simply write to log file if it isn't present
def are_gene_features_present(record, log_h):
    for feature in record.features:
        if feature.type == 'gene':
            return
    log_h.write("No gene annotation features present!!!\n")
    return

# Remove any features classified as tRNA or rRNA from the feature dictionary
def remove_trna_and_rrna_features(record, log_h):
    trna_count = 0
    rrna_count = 0
    remove_flag = 0
    start_pos = []	# List to store first position of tRNA and rRNA features.  There should be no worry in duplicated entries since these are gene sites
    new_features = []	# List where record.features will be updated in
    for feature in record.features:
        if feature.type == 'tRNA':
            trna_count +=1
            start_pos.append(list(feature.location)[0])	# Push first coordinate of location to start_pos list
        if feature.type == 'rRNA':
            rrna_count +=1
            start_pos.append(list(feature.location)[0])

    # For each feature, iterate through and only add gene and CDS features (leaving out gene/tRNA and gene/rRNA)
    for feature in record.features:
        if feature.type == 'source':
            new_features.append(feature)
            continue	# in case rRNA or tRNA started at the first position, we don't want to remove 'source' by accident
        remove_flag = 0
        for i in start_pos:
            #print str(list(feature.location)[0]) + "\t" + str(i) + "\n";
            if list(feature.location)[0] == i:
                remove_flag = 1
                break
        if remove_flag == 0:
            new_features.append(feature)

    log_h.write("There were " + str(trna_count) + " instances of tRNA features removed and " + str(rrna_count) + " instances of rRNA features removed.\n")
    record.features = new_features
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
                    if feature.location.strand == 1:
                        if list(feature.location)[i] > list(feature.location)[i+1]:	# if prev coordinate is larger than next coord
                            #print feature
                            flag = 1
                            if 'locus_tag' in feature.qualifiers:
                                log_h.write("Feature type " + feature.type + " with locus_tag '" + feature.qualifiers['locus_tag'][0] +
    	                        "' has coordinates that run from the end of the circular DNA back to the beginning.  Deleting feature since this may cause issues later on.\n")
    	                    else:	# could not find locus_tag (i.e. misc features section)
    	                        log_h.write("Feature type " + feature.type + " with coordinates " + str(feature.location) +
    	                        " runs from end of the circular DNA back to the beginning.  Deleting feature since this may cause issues later on.\n")
                            break
                    elif feature.location.strand == -1:	# Handle complementary strands
                        if list(feature.location)[i] < list(feature.location)[i+1]:	# if prev coordinate is smaller than next coord
                            #print feature
                            flag = 1
                            if 'locus_tag' in feature.qualifiers:
                                log_h.write("Feature type " + feature.type + " with locus_tag '" + feature.qualifiers['locus_tag'][0] +
                                "' has coordinates that run from the end of the circular DNA back to the beginning.  Deleting feature since this may cause issues later on.\n")
                            else:
    	                        log_h.write("Feature type " + feature.type + " with coordinates " + str(feature.location) +
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

# Make formats to organism name, such as truncating and word-wrap for overwriting the name later
def format_organism_name(organism, log_h):
    org = organism
    MAX_WIDTH = 80	# Maximum width of a line in Genbank
    HEADER_WIDTH = 12	# Width of the header field in Genbank
    TAG = "  ORGANISM"	# Field name with appropriate spacing in front
    max_len = MAX_WIDTH - HEADER_WIDTH

    #print "ORG_NAME_LEN\t" + str(len(org))

	# Getting much of this code from http://biopython.org/DIST/docs/api/Bio.SeqIO.InsdcIO-pysrc.html
    if len(org) > max_len:
        t_org = org[:(max_len - 4)] + "..."
        #print "T_ORG\t" + t_org
        org_lines = split_name_over_lines(org, max_len)
        n_org = ("\n" + " " * HEADER_WIDTH).join(org_lines)
    	#print "N_ORG\t" + n_org
    	log_h.write("Organism field is too long to write to one line.  Will wrap over to next line")
    else:
    	t_org = org
    	n_org = org

    new = TAG.ljust(HEADER_WIDTH) + n_org + "\n" # Left justify TAG, and add word-wrapped organism line(s) after
    truncated = TAG.ljust(HEADER_WIDTH) + t_org + "\n"	# Left justify TAG, and add truncated organism line after
    #print "NEW\n" + new
    #print "TRUNC\n" + truncated
    return new, truncated

# Takes a given field and splits it into word-wrapped multiple lines given a certain field max width
def split_name_over_lines(field, max_len):
    words = field.split()	# Split field into words to prevent splliting a word across 2 lines
    text = ""
    while words and len(text) + 1 + len(words[0]) <= max_len:	# Use words to build line as close to max_len as possible
        text += " " + words.pop(0)	# Build the line while popping off word list
        text = text.strip()
    answer = [text]	# Assign each formatted line to a list
    while words:	# Are there enough words to fill another line?
        text = words.pop(0)
        while words and len(text) + 1 + len(words[0]) <= max_len:
            text += " " + words.pop(0)
            text = text.strip()
        answer.append(text)
    assert not words	# Sanity check to ensure no words were left behind
    return answer


# Using our up-to-date Genbank record information to write a new Genbank file
def write_output(record_list, outfile):
    out_h = open(outfile, "w")
    SeqIO.write(record_list, out_h, "genbank")
    out_h.close()
    return

# Overwrite fixed organism information to the newly written Genbank file
def overwrite_organism(outfile, organism):

	#  There is 1 pratfall with this method.  If we are looking at a multi-contig Genbank file that
    #  just so happens to have different organsim names for each contig (which shouldn't happen),
    #  and the organism truncation from the output file is the exact same across multiple organism
	#  names, then there is a great chance an incorrect new organism name will be applied

    HEADER_WIDTH = 12
    temp = re.sub("\.gbk", ".tmp", outfile)
    os.rename(outfile, temp)	#move our output file to a temp path
    in_h = open (temp, "r")
    out_h = open(outfile, "w")
    for line in in_h:
        for org in organism:
            org_edits = organism[org]
            if org_edits[1] == line:	# compare truncated line to currently read line
			    line = org_edits[0]	# replace with word_wrapped multiline version of Organism field
	    out_h.write(line)
    in_h.close
    out_h.close

# This will move the final Genbank file to the output directory and remove the dirs and files that
# ultimately will not be used downstream
def delete_extra_content(outdir, base_gbk, log_h):
    preval = outdir + "/prevalidate/"
    preval_gbk = preval + base_gbk + ".gbk"
    val = outdir + "/validate/"
    val_tmp = val + base_gbk + ".tmp"
    val_gbk = val + base_gbk + ".gbk"	# Keeping this, moving up a directory
    final_gbk = outdir + base_gbk + ".gbk"
    os.remove(preval_gbk)
    os.remove(val_tmp)
    os.rename(val_gbk, final_gbk)
    sys.stdout.write("Moved " + val_gbk + " to " + final_gbk + "\n")
    os.rmdir(preval)
    os.rmdir(val)
    sys.stdout.write("Removed unnecessary files and directories\n")
    return

#######
#   MAIN   #
#######

def main():
    # Set up options parser and help usage statement
    usage = "usage: %prog -g /path/to/file.gbk -o /path/to/out/dir"
    description = "Validate a Genbank file. Requires Biopython-1.62 to run"
    parser = OptionParser(usage=usage, description=description)
    parser.add_option("-g", "--genbank_file", help="path to a Genbank-formatted file")
    parser.add_option("-o", "--output_path", help="directory path to write output")
    parser.add_option("-l", "--log_file", help="name of the changelog file (only the name)")
    parser.add_option("-d", "--debug", help="run in debug mode", action="store_true", default=False)
    (options, args) = parser.parse_args()

    if not options.genbank_file:
        parser.error( "Genbank file path (-g) must be provided")
    if not options.output_path:
        parser.error("Output directory path (-o) must be provided")

    gbk = options.genbank_file
    gbk = gbk.rstrip()
#    base_gbk = re.sub("\..+", "", basename(gbk))
    base_gbk = splitext(basename(gbk))[0]
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

    # Initialize logging
    if options.log_file:
        log = outdir + "/" + options.log_file
    else:
        log = outdir + "/" + base_gbk + ".gb_changelog.txt"
    log_h = open(log, "w")
    if options.debug == True:
        sys.stdout.write("DEBUG mode enabled\n")

    pattern = re.compile("\.gb\w*")
    match = pattern.search(gbk)
    if not match:
        sys.stderr.write("File " + gbk + " does not have a proper Genbank file extension (.gbk or .gb)... skipping\n")
        sys.exit(1)
    log_h.write("Now preparing " + gbk + " for validation\n")
    new_gbk = prevalidation(gbk, pre, log_h)
    log_h.write("Now validating " + gbk + " ...\n")
    validate_genbank(new_gbk, val, log_h)

    sys.stdout.write("Finished validating " + gbk + "\n")

    # Leave extra files intact if debug mode is enabled
    if options.debug == False:
        delete_extra_content(outdir, base_gbk, log_h)

    log_h.write("\n---------------\n\n")
    log_h.close()

if __name__ == '__main__':
    main()
    sys.exit(0)

