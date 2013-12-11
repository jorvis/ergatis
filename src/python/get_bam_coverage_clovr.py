#!/usr/bin/env python

"""From a BAM file and reference genome
calculate the amount of genome covered
by a minimum coverage value"""

from __future__ import division
from optparse import OptionParser
from Bio import SeqIO
import subprocess
import sys
import os

class cd:
    """Context manager for changing the current working directory"""
    def __init__(self, newPath):
        self.newPath = newPath

    def __enter__(self):
        self.savedPath = os.getcwd()
        os.chdir(self.newPath)

    def __exit__(self, etype, value, traceback):
        os.chdir(self.savedPath)


def get_seq_length(ref):
    """uses BioPython to get the reference sequence length"""
    infile = open(ref, "rU")
    outfile = open("tmp.txt", "w")
    for record in SeqIO.parse(infile, "fasta"):
	print >> outfile,record.id,len(record.seq)

def get_coverage(bam,size):
    subprocess.check_call("genomeCoverageBed -d -ibam %s -g %s > tmp.out" % (bam,size), shell=True)

def remove_column(temp_file):
    infile = open(temp_file, "rU")
    outfile = open("coverage.out", "w")
    my_fields = [ ]
    for line in infile:
	fields=line.split()
	del fields[1]
	my_fields.append(fields)
    for x in my_fields:
	print >> outfile, "\t".join(x)
	
def sum_coverage(coverage,cov):
    infile = open(coverage, "rU")
    outfile = open("amount_covered.txt", "w")
    all = [ ]
    dict = {}
    for line in infile:
	fields=line.split()
	fields = map(lambda s: s.strip(), fields)
	all.append(fields)
    for x, y in all:
	if int(y)>int(cov):
	   try:
		dict[x].append(y)
	   except KeyError:
		dict[x] = [y]
	else:
	   continue
    for k,v in dict.iteritems():
	print >> outfile, k,"\t",len(v)

def doc(coverage, genome_size, name, suffix):
    incov = open(coverage, "U")
    ingenom = open(genome_size, "U")
    outfile = open("%s_%s_depth.txt" % (name, suffix), "w")
    all = [ ]
    my_dict = {}
    for line in incov:
	fields=line.split()
	fields = map(lambda s: s.strip(), fields)
	all.append(fields)
    for x, y in all:
	if int(y)>int(1):
	   try:
		my_dict[x].append(y)
	   except KeyError:
		my_dict[x] = [y]
	else:
	   continue
    new_dict={}
    for k,v in my_dict.iteritems():
	ints = map(int, v)
	new_dict.update({k:sum(ints)})
    genome_size_dict = {}
    for line in ingenom:
	fields = line.split()
	genome_size_dict.update({fields[0]:fields[1]})
    print >> outfile, name,"\n",
    for k,v in new_dict.iteritems():
	print >> outfile, k,"\t",round(int(v)/int(genome_size_dict.get(k)),0)
    for y,z in genome_size_dict.iteritems():
	if y not in new_dict:
		print >> outfile, y,"\t","0"
    
def merge_files_by_column(column, file_1, file_2, out_file):
    """Takes 2 file and merge their columns based on the column. It is assumed
    that line ordering in the files do not match, so we read both files into memory
    and join them"""
    join_map = {}
    for line in open(file_1):
        row = line.split()
        column_value = row.pop(column)
        join_map[column_value] = row

    for line in open(file_2):
        row = line.split()
        column_value = row.pop(column)
        if column_value in join_map:
            join_map[column_value].extend(row)

    fout = open(out_file, 'w')
    for k, v in join_map.iteritems():
        fout.write('\t'.join([k] + v) + '\n')		

def get_seq_name(in_fasta):
    """used for renaming the sequences"""
    return os.path.basename(in_fasta)

def report_stats(results, bam, name, suffix):
    infile = open(results, "rU")
    outfile = open("%s_%s_breadth.txt" % (name, suffix), "w")
    print >> outfile, name,"\n",
    for line in infile:
	fields = line.split()
	chromosome = fields[0]
	try:
	    amount = (int(fields[2])/int(fields[1]))*100
	    print >> outfile,chromosome,"\t",amount,"\n",
	except:
	    print >> outfile, chromosome,"\t","0","\n",
	    sys.exc_clear()
    infile.close()
    outfile.close()

def test_file(option, opt_str, value, parser):
    try:
        with open(value): setattr(parser.values, option.dest, value)
    except IOError:
        print '%s file cannot be opened' % option
        sys.exit()

def main(bam, ref, cov, out_dir, name, suffix):
    with cd(out_dir):
        get_seq_length(ref)
        subprocess.check_call('sed "s/ /\\t/g" tmp.txt > genome_size.txt', shell=True)
        get_coverage(bam,"genome_size.txt")
        remove_column("tmp.out")
        sum_coverage("coverage.out",cov)
        merge_files_by_column(0,"genome_size.txt", "amount_covered.txt", "results.txt")
        report_stats("results.txt", bam, name, suffix)
        doc("coverage.out", "genome_size.txt", name, suffix)
        subprocess.check_call("rm tmp.txt tmp.out amount_covered.txt coverage.out", shell=True)
	
if __name__ == "__main__":
    usage="usage: %prog [options]"
    parser = OptionParser(usage=usage) 
    parser.add_option("-b", "--input_bam", dest="bam",
                      help="/path/to/bam [REQUIRED]",
                      action="callback", callback=test_file, type="string")
    parser.add_option("-r", "--reference", dest="ref",
                      help="/path/to/reference fasta [REQUIRED]",
                      action="callback", callback=test_file, type="string")
    parser.add_option("-c", "--coverage", dest="cov",
                      help="minimum coverage across the genome",
                      default="10", type="int")
    parser.add_option("-n", "--name", dest="name",
                      help="output prefix name [REQUIRED]",
		      action="store", type="string")
    parser.add_option("-s", "--suffix", dest="suffix",
                      help="output suffix name [REQUIRED]",
		     action="store", type="string")
    parser.add_option("-o", "--output_dir", dest="out_dir", 
                      help="Directory to write files too", type="string")
    options, args = parser.parse_args()
    
    mandatories = ["bam", "ref", "name", "suffix"]
    for m in mandatories:
	if not options.__dict__[m]:
		print "\nMust provide %s.\n" %m
		parser.print_help()
		exit(-1)

    main(options.bam, options.ref, options.cov, options.out_dir, options.name, options.suffix)

