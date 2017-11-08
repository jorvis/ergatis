#!/usr/bin/env julia

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

# Determine number of comparisons to use in pangenome bootstrapping
function determine_comparisons(comp::Integer, mult::Integer, hp::HitsProfile)
	if comp > 0
		return comp
	else
		if mult > 0
			return estimate_comparisons(hp.num_genomes, mult)
		else
			return hp.num_genomes * 1000
		end
	end
end

# Estimate the number of comparisons necessary for good bootstrapping
function estimate_comparisons(genomes::Integer, mult::Integer)
    total_comp = 0
    for i in 2:(genomes+1)
        theor = get_theoretical_comps(genomes, i)
        real = mult * genomes
		total_comp += min(theor, real)
	end
    return total_comp
end

# Get number of combinations for the given size pangenome
function get_theoretical_comps(genomes::Integer, i::Integer)
    factorial(genomes) / (factorial(genomes - i) * factorial(i-1))
end

#############
# Pangenome #
#############

# Create the pangenome database using the HitsProfile matris
function create_dataset(pg::Pangenome, hp::HitsProfile, pangenome_genes_flag)
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

# Get number of permutations for the given size pangenome
function get_permutations(genomes::Integer, i::Integer)
    factorial(genomes) / factorial(genomes - i)
end

# Get number of combinations for the given size pangenome
function get_combinations(genomes::Integer, i::Integer)
	return factorial(genomes) / (factorial(genomes - i) * factorial(i))

# Creates the pangenome subset of a passed in size
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
				write(gene_fh, "---\n")
			end
		end
	end
end

# Will calculate the pangenome size for a particular set of genomes
function calculate_pangenome(genomes::Integer, hp::HitsProfile, gene_fh, print_genes_flag)
	done_genomes = []
	pangenome_size = 0

	for g in genomes
		for gene in keys(hp.profile[g])
			shared = false
			# Increment pangenome size if gene hit isn't shared with a processed genome
			for g2 in done_genomes:
				if haskey(hp.profile[g][gene], g2):
					shared = true
					break
				end
			end
			if not shared:
				pangenome_size += 1
				if print_genes_flag
					write(gene_fh,"$g:$gene\n")
				end
			end
		append!(done_genomes,g)
	end
	return pangenome_size
end

# Write to output file
function write_output(pan_size, set_size, genome_string, fh)
	write(fh, "$set_size\t$pan_size\t$genome_string\n")
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
			default = 0
			range_tester = (x -> x>=0)
		"--multiplicity", "-m"
			help = "Another option for sampling based on a multiplicity factor ((sum(m*n) for n=[2..n]) number of comparisons)"
			metavar = "20"
			arg_type = Int
			default = 0
			range_tester = (x -> x>=0)
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

	args["comparisons"] = determine_comparisons(args["comparisons"], args["multiplicity"], hp)
	create_dataset( Pangenome(output_file, genes_file, args["respect_order"], args["comparisons"]), hp, args["print_pangenome_genes"] )
end

main()

#=

class Pangenome:

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

=#
