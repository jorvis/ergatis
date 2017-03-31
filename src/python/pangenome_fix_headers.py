#!/usr/bin/env python

"""
pangenome_fix_headers.py - Format pangenome headers so that genomes and genes are easily identifiable and parseable
By: Shaun Adkins (sadkins@som.umarylane.edu)

"""

# vim: tabstop=8 expandtab shiftwidth=4 softtabstop=4

from argparse import ArgumentParser
import sys
from Bio import SeqIO

# Read in a fasta file as a list of SeqRecords
def read_in_fasta(fasta_in):
	#print "reading in fasta file\n"
	return list( SeqIO.parse(fasta_in, "fasta") )

# Parse mapping file and keep specific fields
def parse_map_file(mapping_in):
	#print "parsing mapping file\n"
	lines = []
	with open (mapping_in, 'r') as mapping_fh:
		for line in mapping_fh:
			vals = line.split("\t")
			seq_id, gene = vals[0].split('|||')
# Keeping sequence id, polypeptid id, gene name, and the organism strain name
			kept_vals = (seq_id, vals[5], gene, vals[7])
			lines.append(kept_vals)
	return lines

# Replace the current fasta header with a properly named one
def replace_fasta_headers(records, replace_vals):
	#print "replacing fasta headers\n"
	updated_records = []
	for record in records:
		header = record.id
		#print record.id + "\n";
# Split on first space
# NOTE: This is under the assumption that the bsml2fasta component for assemblies has USE_SEQUENCE_IDS_IN_FASTA set to 1 in the config file.  The polypeptide version does not need this
		fields = header.split(' ', 1)
		seq_id = fields[0]
# Search for seq_id amongst our kept mapping values
		for line in replace_vals:
# Handle the 'assembly' cases or the 'polypeptide' cases
			if seq_id == line[0]:	# Assembly
				record.id = 'gnl|' + line[3] + '|' + line[0]	# gnl|org_strain|genome/contig
				record.description = ''
				updated_records.append(record)
				break
			elif seq_id == line[1]:	# Polypeptide
				record.id = 'gnl|' + line[3] + '|' + line[2]	# gnl|org_strain|locus_tag
				record.description = ''
				updated_records.append(record)
				break
	return updated_records

def write_new_fasta(record_list, fasta_out):
	SeqIO.write(record_list, fasta_out, "fasta")
	return

########
# Main #
########

def main():
# Set up options parser and help statement
	description = "Format pangenome headers so that genomes and genes are easily identifiable and parseable"
	parser = ArgumentParser(description=description)
	parser.add_argument("--mapping_file", "-m", help="Path to the mapping file created by mugsyprep.pl", metavar="/path/to/mugsymap_complete.txt", required=True)
	parser.add_argument("--input_file", "-f", help="Path to a fasta file.  File should have information that can make it easy to map the genome and gene using the mugsymap file", metavar="/path/to/genome.fa", required=True)
	parser.add_argument("--output_file", "-o", help="Path to write the output file", metavar="/path/to/output.txt", required=True)
	args = parser.parse_args()

	replace_vals = parse_map_file(args.mapping_file)
	records = read_in_fasta(args.input_file)
	new_records = replace_fasta_headers(records, replace_vals)
	write_new_fasta(new_records, args.output_file)

if __name__ == '__main__':
	main()
	sys.exit(0)
