#!/usr/bin/julia

# pangenome_make_table.jl - Reads a pangenome profile and generates core, new, and shared gene hits for the pangenome
# By: Shaun Adkins (sadkins@som.umaryland.edu)

using ArgParse

#############
# Type defs #
#############

type HitsProfile
    profile::Dict
    genome_names::Vector{String}
    num_genomes::Integer
end

HitsProfile() = HitsProfile(Dict(), Vector{String}(), 0)

###############
# HitsProfile #
###############

function read_in_profile(hp::HitsProfile, profile)
    # Read in tab-delimited file, and grab genome names
    (data, headers) = readdlm(profile, '\t', header=true, comments=false)
    hp.num_genomes = length(headers) - 2
    hp.genome_names = headers[3:end]
    hp.profile = parse_data(hp.genome_names, data)
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

function determine_sampling_rate( hp::HitsProfile, args)
    if args["comparisons"] > 0
        return estimate_multiplicity(hp.num_genomes, args["comparisons"])
    else
        return args["multiplicity"]
    end
end

function estimate_multiplicity( genomes::Integer, comparisons::Integer)
    # Determine multiplex sampling rate based on specified comparisons
    lower_mult = 0
    lower_comp = 0
    lower_theor = 0
    upper_mult = 0
    upper_comp = 0
    upper_theor = 0
    ldiff = 0
    udiff = 0

	# Bound limits for estimating multiplicity rate
	lower_bounds = 5
	upper_bounds = 5000

    # Establish multiplex bounds
    for i = lower_bounds:upper_bounds
        (est_comps, theor_comps) = estimate_comparisons(genomes, i)
        if est_comps < comparisons
            lower_mult = i
            lower_comp = est_comps
            lower_theor = theor_comps
        else
            upper_mult = i
            upper_comp = est_comps
            upper_theor = theor_comps
            break
        end
    end
    if upper_mult == lower_bounds
        (est_comps, theor_comps) = estimate_comparisons(genomes, lower_bounds)
    elseif lower_mult == upper_bounds
        return upper_bounds
    else
        lower_diff = comparisons - lower_comp
        upper_diff = upper_comp - comparisons
        if ldiff <= udiff
            return lower_mult
        else
            return upper_mult
        end
    end
end

function estimate_comparisons(genomes::Integer, mult::Integer)
    # Estimate the number of comparisons necessary for good bootstrapping
    total_comp = 0
    theor_comp = 0
    for i = range(2,genomes+1)
        theor = get_theoretical_comps(genomes, i)
        real = mult * genomes
        theor_comp += theor
        total_comp += min(theor, real)
    end
    return (total_comp, theor_comp)
end

function get_theoretical_comps(genomes::Integer, i::Integer)
    # Get number of combinations for the given size pangenome
    factorial(genomes) / (factorial(genomes - i) * factorial(i-1))
end

function create_reference_string(reference_set)
    # Create a string from a random set of genomes
    seen_vec = [ reference_set[i] for i = 1:length(reference_set) ]
    sort!(seen_vec)
    join(seen_vec, ":")
end

function categorize_shared_gene(genes_by_category, comp_genome::String, gene::String)
    # Increment shared count for this gene if we have a hit
    get!(genes_by_category[comp_genome]["shared"], gene, 0)
    genes_by_category[comp_genome]["shared"][gene] += 1
    return genes_by_category
end

function categorize_core_gene(genes_by_category, comp_genome::String, gene::String)
    # Determine if the given gene is present in all other genomes up to this point.
    genes_by_category[comp_genome]["core"][gene] = 1
    return genes_by_category
end

function categorize_new_gene(genes_by_category, comp_genome::String, gene::String)
    # Determine if the given gene is present in no other genomes up to this point.
    genes_by_category[comp_genome]["new"][gene] = 1
    return genes_by_category
end

function write_single_genomes(hp::HitsProfile, fh)
    # Write the number of genes for the individual genomes
    genome_index = Dict{String, Integer}(hp.genome_names[i] => i for i in 1:hp.num_genomes)

    # For each genome, get number of genes and write to file (also write assigned index)
    for genome in hp.genome_names
        num_genes = length(keys(hp.profile[genome]))
        index = genome_index[genome]
        write(fh, "1\t0\t0\t$num_genes\t$index\t$genome\n")
    end
end

function write_results(genes_by_category, comp_genome::String, size::Integer, hp::HitsProfile, fh)
    index = findfirst(hp.genome_names,comp_genome)
    p_size = size + 1 # Pangenome size = Ref genomes + seed genome

    # Write core/shared/new results to file
    shared_count = length(keys(genes_by_category[comp_genome]["shared"]))
    new_count = length(keys(genes_by_category[comp_genome]["new"]))
    core_count = length(keys(genes_by_category[comp_genome]["core"]))

    write(fh, "$p_size\t$core_count\t$shared_count\t$new_count\t$index\t$comp_genome\n")
end

function do_analysis_with_sampling(args, hp::HitsProfile, fh)
    multiplex = determine_sampling_rate(hp, args)

    # Iterate through each possible pangenome size
    for i = 1:(hp.num_genomes-1)
        tic()
        print("Running analysis with $i genomes...")
        sampling_counter = Dict{String, Integer}()

        fac_numer = factorial(big(hp.num_genomes-1))    # Must convert to BigInt first
        fac_denom1 = factorial(big((i+1)-1))            # ... otherwise OverflowError
        fac_denom2 = factorial(big(hp.num_genomes-(i+1)))
        true_max = fac_numer / (fac_denom1 * fac_denom2)

        # Iterate through the all genomes in set
        for j = 1:hp.num_genomes
            seen_strings = Dict{String, Bool}()
            seed_genome = hp.genome_names[j]
            # Create group of reference genomes minus the seed genome
            ref_genomes = copy(hp.genome_names) # copying the variable creates a pointer
            splice!(ref_genomes, j)
            point_count = 0
            while point_count < min(multiplex, true_max)
                # Create of string of N-reference genomes for an N+1 sized pangenome
                reference_set = rand(ref_genomes, i)
                sort_string = create_reference_string(reference_set)

                # Proceed if this is a new string
                if haskey(seen_strings, sort_string) == false
                    get!(sampling_counter, seed_genome, 1)
                    genes_by_category = Dict()
                    genes_by_category[seed_genome] = Dict("shared" => Dict(), "new" => Dict(), "core" => Dict())

                    # Iterate through all the genes in the seed genome
                    for gene in keys(hp.profile[seed_genome])
                        for ref_g in reference_set
                            # If reference genome is a key in the HitsProfile, add to "shared" genes count
                            if haskey(hp.profile[seed_genome][gene], ref_g)
                                genes_by_category = categorize_shared_gene(genes_by_category, seed_genome, gene)
                            end
                        end
                        # Determine if gene is shared by all or none of the referece genomes
                        count = get(genes_by_category[seed_genome]["shared"], gene, 0)
                        if count == length(reference_set)
                            genes_by_category = categorize_core_gene(genes_by_category, seed_genome, gene)
                        elseif count == 0
                            genes_by_category = categorize_new_gene(genes_by_category, seed_genome, gene)
                        end
                    end
                    write_results(genes_by_category, seed_genome, i, hp, fh)
                    sampling_counter[seed_genome] += 1
                    seen_strings[sort_string] = true
                    point_count += 1
                end
            end
        end

        elapsed = toq()
        println("time elapsed: $elapsed secs")
    end

end

function do_analysis_without_sampling()
    println("Doing analysis without sampling is not implemented at this time")
end

function parse_commandline()
    s = ArgParseSettings(description = "Reads a pangenome profile matrix and generates total genes in the pangenome",
    prog = "pangenome_make_pangenome.jl",
	version = "1.0",
	add_version = true,
	add_help = true)

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
			default = 0
			range_tester = (x->x>=0)
            arg_type = Int
        "--multiplicity", "-m"
            help = "Another option for sampling based on a multiplicity factor ((sum(m*n) for n=[2..n]) number of comparisons)"
            metavar = "20"
			default = 0
			range_tester = (x->x>=0)
            arg_type = Int
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
    output_file = "$output_dir/core_shared_new.tsv"
    out_fh = open(output_file, "w")

    ## our output table will take the form:
    ## genome_count  core    shared    new  root_genome
    write_single_genomes(hp, out_fh)

    if args["comparisons"] == 0 && args["multiplicity"] == 0
        do_analysis_without_sampling()
    else
        do_analysis_with_sampling(args, hp, out_fh)
    end

    close(out_fh)
end

main()
quit()
