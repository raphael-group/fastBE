params.proj_dir      = "/n/fs/ragr-research/projects/fastBE"
params.outputDir     = "/n/fs/ragr-research/projects/fastBE/nextflow_results/"
params.sim_script    = "${params.proj_dir}/scripts/simulation.py"

params.allele_minima = "${params.proj_dir}/build/src/fastbe"

params.nmutations   = [3000]
params.nclones      = [25, 50, 100, 200, 500, 1000]
params.nsamples     = [200]
params.temperatures = [0.01, 0.05, 0.1, 0.25, 0.5]
params.seeds        = [0, 1, 2, 3]
params.coverage     = [100]

process create_sim {
    cpus 1
    memory '4 GB'
    time '59m'
    errorStrategy 'ignore'

    input:
        tuple val(mutations), val(samples), val(clones), val(coverage), val(seed), val(temperature)

    output:
        tuple file("sim_clonal_matrix.txt"), file("sim_mutation_to_clone_mapping.txt"), file("sim_obs_frequency_matrix.txt"), 
              file("sim_total_matrix.txt"), file("sim_tree.txt"), file("sim_usage_matrix.txt"), file("sim_variant_matrix.txt"),
              val(clones), val(temperature), val("m${mutations}_n${clones}_s${samples}_c${coverage}_r${seed}_t${temperature}")

    """
    python '${params.sim_script}' --mutations ${mutations} --samples ${samples} --clones ${clones} --coverage ${coverage} --seed ${seed} --output sim
    """
}

process allele_minima {
    cpus 16
    memory '4 GB'
    time '24h'
    //errorStrategy 'ignore'

    input:
        tuple path(clonal_matrix), path(mut_clone_mapping), path(freq_matrix), path(total_matrix),
              path(clone_tree), path(usage_matrix), path(variant_matrix), val(clones), val(temperature), val(id)

    output:
        tuple file("inferred_tree.txt"), file("inferred_results.json"), file("timing.txt"), file("inferred_sample_0_sim_annealing_info.json"), val(id)

    """
    /usr/bin/time -v '${params.allele_minima}' search ${freq_matrix} --initial_probability ${temperature} --algorithm simulated_annealing \
                  -s 1 --output inferred -t ${task.cpus} --max_iterations 500000 --logging 2>> timing.txt
    """
}

workflow {
    parameter_channel = channel.fromList(params.nmutations)
                               .combine(channel.fromList(params.nsamples))
                               .combine(channel.fromList(params.nclones))
                               .combine(channel.fromList(params.coverage))
                               .combine(channel.fromList(params.seeds))
                               .combine(channel.fromList(params.temperatures))

    parameter_channel = parameter_channel.filter {
      it[0] >= 2*it[2] // require twice as many mutations as clones
    }

    // create directories
    file("${params.outputDir}/simulated_annealing/allele_minima/").mkdirs()

    // run simulation
    simulation = parameter_channel | create_sim

    // save ground truth
    simulation | map { 
        clonal_matrix, mut_clone_mapping, freq_matrix, total_matrix,
        clone_tree, usage_matrix, variant_matrix, clones, temperature, id ->
        
        outputPrefix = "${params.outputDir}/ground_truth/${id}"
        clonal_matrix.copyTo("${outputPrefix}_clonal_matrix.txt")
        mut_clone_mapping.copyTo("${outputPrefix}_mutation_clone_mapping.txt")
        freq_matrix.copyTo("${outputPrefix}_frequency_matrix.txt")
        total_matrix.copyTo("${outputPrefix}_total_matrix.txt")
        variant_matrix.copyTo("${outputPrefix}_variant_matrix.txt")
        usage_matrix.copyTo("${outputPrefix}_usage_matrix.txt")
        clone_tree.copyTo("${outputPrefix}_tree.txt")
    }

    // run AlleleMinima
    simulation | allele_minima | map { inferred_tree, inferred_results, timing, sim_annealing, id ->
        outputPrefix = "${params.outputDir}/simulated_annealing/allele_minima/${id}"
        inferred_tree.moveTo("${outputPrefix}_inferred_tree.txt")
        inferred_results.moveTo("${outputPrefix}_inferred_results.json")
        sim_annealing.moveTo("${outputPrefix}_sim_annealing_results.json")
        timing.moveTo("${outputPrefix}_timing.txt")
    }
}
