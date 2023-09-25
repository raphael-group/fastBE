params.proj_dir      = "/n/fs/ragr-research/projects/vafpp"
params.outputDir     = "/n/fs/ragr-research/projects/vafpp/nextflow_results/benchmarking/"
params.sim_script    = "${params.proj_dir}/scripts/simulation.py"

params.allele_minima = "${params.proj_dir}/build/src/vafpp"

params.nmutations = [200]
params.nclones    = [3, 5, 10, 30, 50]
params.nsamples   = [5, 10, 25, 50, 100]
params.seeds      = [0, 1]
params.coverage   = [100]

process create_sim {
    cpus 1
    memory '4 GB'
    time '59m'
    errorStrategy 'ignore'

    input:
        tuple val(mutations), val(samples), val(clones), val(coverage), val(seed)

    output:
        tuple file("sim_clonal_matrix.txt"), file("sim_mutation_to_clone_mapping.txt"), file("sim_obs_frequency_matrix.txt"), 
              file("sim_total_matrix.txt"), file("sim_tree.txt"), file("sim_usage_matrix.txt"), file("sim_variant_matrix.txt"),
              val(clones), val("m${mutations}_n${clones}_s${samples}_c${coverage}_r${seed}")

    """
    python '${params.sim_script}' --mutations ${mutations} --samples ${samples} --clones ${clones} --coverage ${coverage} --seed ${seed} --output sim
    """
}

process allele_minima {
    cpus 16
    memory '2 GB'
    time '24h'
    errorStrategy 'ignore'

    input:
        tuple path(clonal_matrix), path(mut_clone_mapping), path(freq_matrix), path(total_matrix),
              path(clone_tree), path(usage_matrix), path(variant_matrix), val(clones), val(id)

    output:
        tuple file("inferred_tree.txt"), file("inferred_results.json"), file("timing.txt"), val(id)

    """
    /usr/bin/time -v '${params.allele_minima}' search ${freq_matrix} -a 1 -s 128 --output inferred -t ${task.cpus} 2>> timing.txt
    """
}

workflow {
    parameter_channel = channel.fromList(params.nmutations)
                               .combine(channel.fromList(params.nsamples))
                               .combine(channel.fromList(params.nclones))
                               .combine(channel.fromList(params.coverage))
                               .combine(channel.fromList(params.seeds))

    parameter_channel = parameter_channel.filter {
      it[0] >= 2*it[2] // require twice as many mutations as clones
    }

    run_id = 5
    simulation = parameter_channel | create_sim

    // run AlleleMinima
    simulation | allele_minima | map { inferred_tree, inferred_results, timing, id ->
        outputPrefix = "${params.outputDir}/allele_minima/r${run_id}_${id}"
        inferred_tree.moveTo("${outputPrefix}_inferred_tree.txt")
        inferred_results.moveTo("${outputPrefix}_inferred_results.json")
        timing.moveTo("${outputPrefix}_timing.txt")
    }
}
