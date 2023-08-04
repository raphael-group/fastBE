params.proj_dir      = "/n/fs/ragr-research/projects/vafpp"
params.sim_script    = "${params.proj_dir}/scripts/simulation.py"
params.python_script = "${params.proj_dir}/scripts/vafpp_lp.py"
params.cpp_command   = "${params.proj_dir}/build/src/vafpp"
params.outputDir     = "/n/fs/ragr-research/projects/vafpp/nextflow_results/"

params.nmutations = [20, 60, 100, 200]
params.nclones    = [10, 30, 50, 100]
params.nsamples   = [3, 10, 25, 50, 100]
params.seeds      = [0, 1, 2, 3, 4, 5]
params.coverage   = [100]

process create_sim {
    cpus 1
    memory '4 GB'
    time '59m'

    input:
        tuple val(mutations), val(samples), val(clones), val(coverage), val(seed)

    output:
        tuple file("sim_clonal_matrix.txt"), file("sim_mutation_to_clone_mapping.txt"), file("sim_obs_frequency_matrix.txt"), 
              file("sim_total_matrix.txt"), file("sim_tree.txt"), file("sim_usage_matrix.txt"), file("sim_variant_matrix.txt"),
              val("m${mutations}_n${clones}_s${samples}_c${coverage}_r${seed}")

    """
    python '${params.sim_script}' --mutations ${mutations} --samples ${samples} --clones ${clones} --coverage ${coverage} --seed ${seed} --output sim
    """
}

process infer_cpp {
    cpus 16
    memory '2 GB'
    time '24h'

    input:
        tuple path(clonal_matrix), path(mut_clone_mapping), path(freq_matrix), path(total_matrix),
              path(clone_tree), path(usage_matrix), path(variant_matrix), val(id)

    output:
        tuple file("inferred_tree.txt"), file("inferred_results.json"), val(id)

    """
    '${params.cpp_command}' search ${freq_matrix} -a 0.15 -s 200 --output inferred -t ${task.cpus}
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

    create_sim([50, 10, 10, 100, 0]) | infer_cpp | map { inferred_tree, inferred_results, id ->
        outputPrefix = "${params.outputDir}/cpp/${id}"
        inferred_tree.moveTo("${outputPrefix}_inferred_tree.txt")
        inferred_results.moveTo("${outputPrefix}_inferred_results.json")
    }
}
