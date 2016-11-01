#!/usr/bin/python

"""
pangenome_make_pangenome.py - Reads a pangenome profile and generates total genes in the pangenome
By: Shaun Adkins (sadkins@som.umaryland.edu)

"""

from argparse import ArgumentParser
import sys
import logging
import re
from math import factorial
from math import floor
import random

#############
# Pangenome #
#############

class Pangenome:

    def __init__(self, out, ro, comp, mult, hp):
        self.logger = logging.getLogger(__name__)
        self.logger.info("Creating a Pangenome instance")
        self.output_file = out + "/pangenome.output"
        self.genes_file = self.output_file + ".genes"
        self.respect_order = ro
        self.comparisons = comp if comp > 0 else hp.num_genomes * 1000
        if mult > 0:
            est_comp = self.estimate_comparisons(hp.num_genomes, mult)
            self.comparisons = est_comp
        self.logger.info("Estimated number of comparisons is " + str(self.comparisons))


    def estimate_comparisons(self, genomes, mult):
        """ Estimate the number of comparisons necessary to make a good pangenome """
        self.logger.info("Estimating the number of comparisons")
        total_comp = 0
        for i in range(2,genomes+1):
            theor = self.get_theoretical_comps(genomes, i)
            real = mult * genomes
            if theor < real:
                total_comp += theor
            else:
                total_comp += real
        return total_comp

    def get_theoretical_comps(self, genomes, i):
        """ Calculate all possible comparisons """
        return factorial(genomes) / (factorial(genomes - i) * factorial(i-1))

    def create_dataset(self, hp, pangenome_genes_flag):
        """ Create the pangenome database using the profile matrix """
        self.logger.info("Beginning to create the pangenome")
        self.logger.info("Respecting dataset order: " + str(self.respect_order))
        out_fh = self.open_for_writing(self.output_file)
        gene_fh = self.open_for_writing(self.genes_file)
        for i in range(1, hp.num_genomes+1):
            self.logger.info("Now calulating pangenome for size " + str(i))
            if self.respect_order:
                max1 = floor(self.get_permutations(hp.num_genomes, i) + 0.5)
            else:
                max1 = floor(self.get_combinations(hp.num_genomes, i) + 0.5)
            max2 = floor((self.comparisons/hp.num_genomes) + 0.5)
            # Only print pangenome genes for maximum num of genomes
            if pangenome_genes_flag and i == hp.num_genomes:
                print_genes_flag = 1
            else:
                print_genes_flag = 0
            self.select_genomes(max1, max2, hp, i, out_fh, gene_fh, print_genes_flag)
        gene_fh.close()
        out_fh.close()
        self.logger.info("Finished writing the pangenome output")

    def get_permutations(self, genomes, i):
        """ Get number of permutations for the given size pangenome """
        return factorial(genomes) / factorial(genomes - i)

    def get_combinations(self, genomes, i):
        """ Get number of combinations for the given size pangenome """
        return factorial(genomes) / (factorial(genomes - i) * factorial(i))

    def select_genomes(self, max1, max2, hp, size, out_fh, gene_fh, print_genes_flag):
        """ Creates the pangenome subset of a passed in size """
        iteration = 0
        seen = []
        while iteration < max1 and iteration < max2:
            genome_set = random.sample(hp.genome_names, size)
            if self.respect_order:
                random.shuffle(genome_set)
            genome_string = '-'.join(genome_set)
            if genome_string not in seen:
                seen.append(genome_string)
                pangenome_size = self.calculate_pangenome(genome_set, hp, gene_fh, print_genes_flag)
                self.write_output(pangenome_size, size, genome_string, out_fh)
                iteration += 1

                if print_genes_flag:
                    gene_fh.write("---\n")

    def calculate_pangenome(self, genomes, hp, gene_fh, print_genes_flag):
        """ Will calculate the pangenome size for a particular set of genomes """
        done_genomes = []
        pangenome_size = 0

        for g in genomes:
            for gene in hp.profile[g].keys():
                shared = False
                # Increment pangenome size if gene hit isn't shared with a processed genome
                ### I tried to implement any(x in list for x in dict) but it was slower
                for g2 in done_genomes:
                    if g2 in hp.profile[g][gene]:
                        shared = True
                        break
                if not shared:
                    pangenome_size += 1
                    if print_genes_flag:
                        gene_fh.write(g + ":" + gene + "\n")
                    self.logger.debug(gene + " from genome " + g + " was added to the pangenome count")
            done_genomes.append(g)
        return pangenome_size

    def write_output(self, pan_size, set_size, genome_string, fh):
        """ Write to output file """
        fh.write(str(set_size) + "\t" + str(pan_size) + "\t" + genome_string + "\n")

    def open_for_writing(self, file):
        """ Open a file for writing """
        try:
            fh = open(file, "w")
        except IOError:
            self.logger.error("Cannot open file " + file + " for writing!")
            sys.exit(1)
        else:
            return fh

##########
# HitsProfile #
##########

class HitsProfile:

    def __init__(self):
        self.logger = logging.getLogger(__name__)
        self.logger.info("Creating a HitsProfile instance")
        self.genome_names = []
        self.num_genomes = 0
        self.profile = {}

    def read_in_profile(self, profile):
        """ Read in the profile hit matrix file """
        self.logger.info("Preparing to read in " + profile)
        profile_fh = self.open_file(profile)
        self.read_first_line(profile_fh)
        self.read_in_hits(profile_fh)
        profile_fh.close()

    def open_file(self, file):
        """ Open a file for reading """
        try:
            fh = open(file, "r")
        except IOError:
            self.logger.error("Cannot open file " + file + " for reading!")
            sys.exit(1)
        else:
            return fh

    def read_first_line(self, fh):
        """ Parse first line of the profile file """
        self.logger.info("Parsing first line of profile")
        line = fh.readline()
        line = line[2:] # Chop off first two characters (#\s)
        gen_names = re.split(r'\t+', line.rstrip())
        self.num_genomes =len(gen_names) - 2
        self.genome_names = self.index_genome_names(gen_names)
        self.logger.debug(self.genome_names)

    def index_genome_names(self, gen_names):
        """ Creates a list of tuples specifying the index number and the genome name """
        self.logger.info("Indexing individual genome names")
        # Genome names start at i=2, but the new list will be at i=0
        return [ gen_names[i] for i in range(2, len(gen_names)) ]

    def read_in_hits(self, fh):
        """ Reads in the hits from the pangenome profile matrix """
        self.logger.info("Storing gene hits per genome")
        for line in fh:
            gene_line = re.split(r'\t+', line.rstrip('\t'), 2)
            genome, gene, hits = gene_line
            h = re.split(r'\t+', hits)
            #self.logger.debug(hits)
            genome_hits = self.filter_gene_hits(h)
            # If genome isn't in dict, initialize and add gene/hits
            try:
                self.profile[genome][gene] = genome_hits
            except KeyError:
                self.logger.debug("Initializing genome " + genome + " into hash")
                self.profile[genome] = {gene : genome_hits}

            self.logger.debug(genome + " - " + gene + " - " + str(len(genome_hits)))

    def filter_gene_hits(self, hits):
        """ Get all genomes with a hit in the current gene """
        return {self.genome_names[i] : 1 for i in range(len(hits)) if int(hits[i]) == 1}

########
# Main #
########

def main():
    # Set up options parser and help statement
    description = "Reads a pangenome profile matrix and generates total genes in the pangenome"
    parser = ArgumentParser(description=description)
    parser.add_argument("--profile", "-p", help="Path to profile matrix", metavar="/path/to/pangenome/profile.txt", required=True)
    parser.add_argument("--output_path", "-o", help="Path to write the output file", metavar="/path/to/pangenome/dir/", required=True)
    parser.add_argument("--comparisons", "-c", help="The number of comparisons to make for any 1 value of N (sampling)", metavar="100000", type=int)
    parser.add_argument("--multiplicity", "-m", help="Another option for sampling based on a multiplicity factor ((sum(m*n) for n=[2..n]) number of comparisons)", metavar="20", type=int)
    parser.add_argument("--respect_order", "-r", help="If enabled, will use permutations instead of combinations", action="store_true", default=False)
    parser.add_argument("--print_pangenome_genes", "-g", help="Enable to output the gene names and genomes for each maximum-size pangenome", action="store_true", default=False)
    parser.add_argument("--log_file", "-l", help="Path to write the logfile", metavar="/path/to/logfile.log")
    parser.add_argument("--debug", "-d", help="Set the debug level", default="ERROR", metavar="DEBUG/INFO/WARNING/ERROR/CRITICAL")
    args = parser.parse_args()
    check_args(args, parser)

	if not args.comparisons:
		args.comparisons = 0
	if not args.multiplicity:
		args.multiplicity = 0

    # Instantiate a HitsProfile class
    hp = HitsProfile()
    hp.read_in_profile(args.profile)

    pg = Pangenome(args.output_path, args.respect_order, args.comparisons, args.multiplicity, hp)
    pg.create_dataset(hp, args.print_pangenome_genes)


def check_args(args, parser):
	""" Validate the passed arguments """
	logger = configure_logger(args.log_file, args.debug)
	if not (args.comparisons or args.multiplicity):
		logger.error('Must either specify the "comparisons" or "multiplicity" options')
		sys.exit(1)

def configure_logger(filename, log_level):
    """ Creates a logger object with the appropriate level """
    num_level = getattr(logging, log_level.upper())

	# Verify that our specified log_level has a numerical value associated
    if not isinstance(num_level, int):
        raise ValueError('Invalid log level: %s' % log_level)

	# Create the logger
    logger = logging.getLogger(__name__)
    logger.setLevel(logging.DEBUG)

	# Add console handler
    ch = logging.StreamHandler()
    ch.setLevel(num_level)
    logger.addHandler(ch)

	# If a log_file argument was provided, write to that too
    if filename:
        log_fh = logging.FileHandler(filename, mode='w')
        log_fh.setLevel(logging.DEBUG)	# Let's write all output to the logfile
        formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
        log_fh.setFormatter(formatter)
        logger.addHandler(log_fh)

    return logger

if __name__ == '__main__':
    main()
    sys.exit(0)
