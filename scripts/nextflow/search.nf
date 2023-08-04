params.proj_dir      = "/n/fs/ragr-research/projects/vafpp"
params.sim_script    = "${params.proj_dir}/scripts/simulation.py"
params.pairtree_input_script = "${params.proj_dir}/scripts/make_pairtree_input.py"
params.pairtree_parse_output = "${params.proj_dir}/scripts/parse_pairtree_output.py"
params.allele_minima = "${params.proj_dir}/build/src/vafpp"
params.outputDir     = "/n/fs/ragr-research/projects/vafpp/nextflow_results/"
params.pairtree_bin  = "${params.proj_dir}/dependencies/pairtree/bin/pairtree"

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

process allele_minima {
    cpus 16
    memory '2 GB'
    time '24h'

    input:
        tuple path(clonal_matrix), path(mut_clone_mapping), path(freq_matrix), path(total_matrix),
              path(clone_tree), path(usage_matrix), path(variant_matrix), val(id)

    output:
        tuple file("inferred_tree.txt"), file("inferred_results.json"), val(id)

    """
    '${params.allele_minima}' search ${freq_matrix} -a 0.15 -s 500 --output inferred -t ${task.cpus}
    """
}

process create_pairtree_input { 
    cpus 1
    memory '1 GB'
    time '10m'

    input:
        tuple path(clonal_matrix), path(mut_clone_mapping), path(freq_matrix), path(total_matrix),
              path(clone_tree), path(usage_matrix), path(variant_matrix), val(id)

    output:
        tuple file("pairtree_mutations.ssm"), file("pairtree_params.json"), val(id)

    """
    python '${params.pairtree_input_script}' ${variant_matrix} ${total_matrix} ${mut_clone_mapping} -o pairtree
    """
}

process pairtree {
    cpus 16
    memory '8 GB'
    time '24h'

    input:
        tuple file(ssm), file(params_json), val(id)

    output:
        tuple file("results.npz"), file("best_tree.txt"), val(id)

    """
    '${params.pairtree_bin}' --params ${params_json} ${ssm} results.npz
    python '${params.pairtree_parse_output}' results.npz --output best_tree.txt
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

    // create directories
    file("${params.outputDir}/allele_minima/").mkdirs()
    file("${params.outputDir}/pairtree/").mkdirs()
    file("${params.outputDir}/ground_truth/").mkdirs()

    // run simulation
    simulation = create_sim([50, 10, 10, 100, 0])

    // save ground truth
    simulation | map { 
        clonal_matrix, mut_clone_mapping, freq_matrix, total_matrix,
        clone_tree, usage_matrix, variant_matrix, id ->
        
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
    simulation | allele_minima | map { inferred_tree, inferred_results, id ->
        outputPrefix = "${params.outputDir}/allele_minima/${id}"
        inferred_tree.moveTo("${outputPrefix}_inferred_tree.txt")
        inferred_results.moveTo("${outputPrefix}_inferred_results.json")
    }

    // run Pairtree
    simulation | create_pairtree_input | pairtree | map { results, best_tree, id ->
        outputPrefix = "${params.outputDir}/pairtree/${id}"
        results.moveTo("${outputPrefix}_results.npz")
        best_tree.moveTo("${outputPrefix}_best_tree.txt")
    }
}
