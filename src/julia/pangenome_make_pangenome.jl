#!/usr/bin/julia

#pangenome_make_pangenome.jl - Reads a pangenome profile and generates total genes in the pangenome
#By: Shaun Adkins (sadkins@som.umarylane.edu)

using ArgParse

############
# Type defs #
############

type Pangenome
	output_file::String
	genes_file::String
	respect_order::Bool
	comparisons::Integer
end

type HitsProfile
	profile::Dict
	genome_names::Vector{String}
	num_genomes::Integer
end

HitsProfile() = HitsProfile(Dict(), Vector{String}(), 0)

#############
# Pangenome #
#############

function create_dataset(pg::Pangenome, hp::HitsProfile, pangenome_genes_flag)
	# Create the pangenome database using the HitsProfile matris
	out_fh = open(pg.output_file, "w")
	gene_fh = open(pg.genes_file, "w")
	for i = range(hp.num_genomes)
		if pg.respect_order
			max1 = floor( get_permutations(hp.num_genomes, i) + 0.5)
		else
			max1 = floor( get_combinations(hp.num_genomes, i) + 0.5)
		end
		max2 = floor( pg.comparisons / hp.num_genomes) + 0.5)

		# Only print pangenome genes for max number of genomes if option set
		if pangenome_genes_flag and i == hp.numgenomes
			print_genes_flag = true
		else
			print_genes_flag = false
		end

		select_genomes(max1, max2, hp, i, out_fh, gene_fh, print_genes_flag, pg.respect_order)
	end
	close(out_fh)
	close(gene_fh)
end

function get_permutations(genomes::Integer, i::Integer):
    # Get number of permutations for the given size pangenome
    factorial(genomes) / factorial(genomes - i)
end

function get_theoretical_comps(genomes::Integer, i::Integer):
	# Get number of combinations for the given size pangenome
    factorial(genomes) / (factorial(genomes - i) * factorial(i-1))
end

function estimate_comparisons(genomes::Integer, mult::Integer)
    # Estimate the number of comparisons necessary for good bootstrapping
    total_comp = 0
    for i = range(2,genomes+1)
        theor = get_theoretical_comps(genomes, i)
        real = mult * genomes
		total_comp += min(theor, real)
	end
    return total_comp
end

function select_genomes(max1::Integer, max2::Integer, set_size::Integer, out_fh, gene_fh, print_genes_flag::Bool, respect_order::Bool)
	iteration = 0
	seen = []
	while iterations < max1 and iterations < max2
		genome_set = rand(hp.genome_names, pg_size)
		if respect_order
			shuffle!(genome_set)
		end
		genome_string = join(genome_set, "-")

		if in(genome_string, seen)
			push!(seen, genome_string)
			pangenome_size = calculate_pangenome(genome_set, hp, gene_fh, print_genes_flag)
			write_output(pangenome_size, set_size, genome_string, out_fh)
			iteration += 1

			if print_genes_flag
				#write dashes to gene_fh
			end
		end
	end
end

##########
# HitsProfile #
##########

function read_in_profile(hp::HitsProfile, profile)
	# Read in tab-delimited file, and grab genome names
	(data, headers) = readdlm(profile, '\t', header=true, comments=false)
	hp.num_genomes = length(headers) - 2
	hp.genome_names = headers[3:end]
	hp.profile = parse_data(genome_names, data)
end

function parse_data(headers::Vector{String}, data)
    # Iterate through each line of data
    profile = Dict{String, Any}()
    for i = 1:size(data,1)
        genome = data[i,1]
        gene = data[i,2]
        # Hits to each subject genome
        hits = data[i,3:end]
        g = Dict{String, Bool}()
        for j = 1:length(hits)
            if hits[j] == 1
                g[headers[j]] = true
            end
        end

        # Create typed dictionary for genome if one does not exist
        try
			profile[genome][gene] = g
		catch KeyError e
			profile[genome] = Dict{String, Any}(gene => g)
		end
    end
    return profile
end

########
# Main #
########

function parse_commandline()
    s = ArgParseSettings(description = "Reads a pangenome profile matrix and generates total genes in the pangenome",
						 prog = "pangenome_make_pangenome.jl")

	# The macro to add a table of arguments and options to the ArgParseSettings object
	@add_arg_table s begin
        "--profile", "-p"
            help = "Path to profile matrix"
			metavar = "/path/to/pangenome/profile.txt"
            required = true
		"--output_dir", "-o"
		    help = "Directory to write the output"
			metavar = "/path/to/pangenome/dir/"
			required = true
		"--comparisons", "-c"
			help = "The number of comparisons to make for any 1 value of N (sampling)"
			metavar = "100000"
			arg_type = Int
		"--multiplicity", "-m"
			help = "Another option for sampling based on a multiplicity factor ((sum(m*n) for n=[2..n]) number of comparisons)"
			metavar = "20"
			arg_type = Int
		"--respect_order", "-r"
			help = "If enabled, will use permutations instead of combinations"
			action = :store_true
			default = false
		"--print_pangenome_genes", "-g"
			help = "Enable to output the gene names and genomes for each maximum-size pangenome"
			action = :store_true
			default = false
			# Will not add log_file or debug options for now
    end

	# Converts the ArgParseSettings object into key/value pairs
    return parse_args(s)
end

function main()
    args = parse_commandline()
    println("Parsed args:")
    for (arg,val) in args
        println("  $arg  =>  $val")
    end

	hp = HitsProfile()
	read_in_profile(hp, args["profile"])

	output_dir = args["output_dir"]
	output_file = "$output_dir/pangenome.output"
	genes_file = "$output_file.genes"

    #pg.create_dataset(hp, args.print_pangenome_genes)
	create_dataset( Pangenome(output_file, genes_file, args["respect_order"], args["comparisons"]), hp, args["print_pangenome_genes"] )
end

main()

#=

class Pangenome:

    def __init__(self, out, ro, hp, comp = None, mult = None):
        self.logger = logging.getLogger(__name__)
        self.logger.info("Creating a Pangenome instance")
        self.output_file = out + "/pangenome.output"
        self.genes_file = self.output_file + ".genes"
        self.respect_order = ro
		#Python used in Ergatis (v2.4.3) does not support ternary form, which was previously used
        if comp != None:
            self.comparisons = comp
        else:
            hp.num_genome * 1000
        if mult != None:
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

=#
