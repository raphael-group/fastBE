params.proj_dir      = "/n/fs/ragr-research/projects/fastBE"
params.outputDir     = "/n/fs/ragr-research/projects/fastBE/nextflow_results/"
params.sim_script    = "${params.proj_dir}/scripts/simulation.py"

params.allele_minima = "${params.proj_dir}/build/src/fastBE"

params.pairtree_input_script = "${params.proj_dir}/scripts/make_pairtree_input.py"
params.pairtree_parse_output = "${params.proj_dir}/scripts/parse_pairtree_output.py"
params.pairtree_bin  = "${params.proj_dir}/dependencies/pairtree/bin/pairtree"

params.calder_jar    = "${params.proj_dir}/dependencies/calder/calder.jar"
params.calder_input_script = "${params.proj_dir}/scripts/make_calder_input.py"
params.calder_parse_output = "${params.proj_dir}/scripts/parse_calder_output.py"

params.citup_python = "/n/fs/ragr-data/users/schmidt/miniconda3/envs/citupenv/bin/python"
params.citup_qip = "/n/fs/ragr-data/users/schmidt/miniconda3/envs/citupenv/bin/run_citup_qip.py"
params.citup_input_script = "${params.proj_dir}/scripts/make_citup_input.py"
params.citup_parse_output = "${params.proj_dir}/scripts/parse_citup_output.py"

params.nmutations = [20, 30, 40, 50, 100, 200]
params.nclones    = [3, 5, 10, 30, 50]
params.nsamples   = [5, 10, 25, 50, 100]
params.seeds      = [0, 1, 2, 3, 4, 5]
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
    memory '4 GB'
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

process create_citup_input { 
    cpus 1
    memory '1 GB'
    time '10m'

    input:
        tuple path(clonal_matrix), path(mut_clone_mapping), path(freq_matrix), path(total_matrix),
              path(clone_tree), path(usage_matrix), path(variant_matrix), val(clones), val(id)

    output:
        tuple file("citup_frequency_matrix.txt"), file("citup_clustering.txt"), val(clones), val(id)

    """
    python '${params.citup_input_script}' ${freq_matrix} --output citup_
    """
}

process create_calder_input { 
    cpus 1
    memory '1 GB'
    time '10m'

    input:
        tuple path(clonal_matrix), path(mut_clone_mapping), path(freq_matrix), path(total_matrix),
              path(clone_tree), path(usage_matrix), path(variant_matrix), val(clones), val(id)

    output:
        tuple file("calder_input.tsv"), val(id)

    """
    python '${params.calder_input_script}' ${variant_matrix} ${total_matrix} ${mut_clone_mapping} -o calder_input.tsv
    """
}

process create_pairtree_input { 
    cpus 1
    memory '1 GB'
    time '10m'

    input:
        tuple path(clonal_matrix), path(mut_clone_mapping), path(freq_matrix), path(total_matrix),
              path(clone_tree), path(usage_matrix), path(variant_matrix), val(clones), val(id)

    output:
        tuple file("pairtree_mutations.ssm"), file("pairtree_params.json"), val(id)

    """
    python '${params.pairtree_input_script}' ${variant_matrix} ${total_matrix} ${mut_clone_mapping} -n $clones -o pairtree
    """
}

process pairtree {
    cpus 16
    memory '8 GB'
    time '24h'
    errorStrategy 'ignore'

    input:
        tuple file(ssm), file(params_json), val(id)

    output:
        tuple file("results.npz"), file("best_tree.txt"), file("timing.txt"), val(id)

    """
    /usr/bin/time -v '${params.pairtree_bin}' --parallel 16 --params ${params_json} ${ssm} results.npz 2>> timing.txt
    python '${params.pairtree_parse_output}' results.npz --output best_tree.txt     
    """
}

process calder {
    cpus 16
    memory '8 GB'
    time '24h'
    errorStrategy 'ignore'

    input:
        tuple file(input_tsv), val(id)

    output:
        tuple file("tree.txt"), val(id)

    """
    export LD_LIBRARY_PATH=/n/fs/ragr-research/projects/fastBE/dependencies/calder/glpk-4.65/lib/jni 
    java -jar '${params.calder_jar}' -i ${input_tsv} -o . -v glpk -N
    python '${params.calder_parse_output}' calder_tree1.dot > tree.txt
    """
}

process citup {
    cpus 16
    memory '8 GB'
    time '24h'
    errorStrategy 'ignore'

    input:
        tuple file(freq_matrix), file(clustering), val(clones), val(id)

    output:
        tuple file("tree.txt"), val(id)

    """
    export PATH=/n/fs/ragr-data/users/schmidt/miniconda3/envs/citupenv/bin:$PATH
    '${params.citup_qip}' ${freq_matrix} ${clustering} citup_results.hdf --loglevel INFO --maxjobs 16 --min_nodes ${clones} --max_nodes ${clones} --submit local
    '${params.citup_python}' '${params.citup_parse_output}' citup_results.hdf > tree.txt
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
    file("${params.outputDir}/calder/").mkdirs()
    file("${params.outputDir}/citup/").mkdirs()
    file("${params.outputDir}/ground_truth/").mkdirs()

    // run simulation
    simulation = parameter_channel | create_sim

    // save ground truth
    simulation | map { 
        clonal_matrix, mut_clone_mapping, freq_matrix, total_matrix,
        clone_tree, usage_matrix, variant_matrix, clones, id ->
        
        outputPrefix = "${params.outputDir}/ground_truth/${id}"
        clonal_matrix.copyTo("${outputPrefix}_clonal_matrix.txt")
        mut_clone_mapping.copyTo("${outputPrefix}_mutation_clone_mapping.txt")
        freq_matrix.copyTo("${outputPrefix}_frequency_matrix.txt")
        total_matrix.copyTo("${outputPrefix}_total_matrix.txt")
        variant_matrix.copyTo("${outputPrefix}_variant_matrix.txt")
        usage_matrix.copyTo("${outputPrefix}_usage_matrix.txt")
        clone_tree.copyTo("${outputPrefix}_tree.txt")
    }

    // filter {5 < it[it.size() - 2] && it[it.size() - 2] <= 10}
    // run CITUP
    // simulation  | create_citup_input | citup | map { tree, id ->
        // outputPrefix = "${params.outputDir}/citup/${id}"
        // tree.moveTo("${outputPrefix}_inferred_tree.txt")
    // }

    // run CALDER
    // simulation | create_calder_input | calder | map { tree, id ->
        // outputPrefix = "${params.outputDir}/calder/${id}"
        // tree.moveTo("${outputPrefix}_inferred_tree.txt")
    // }

    // run AlleleMinima
    simulation | allele_minima | map { inferred_tree, inferred_results, timing, id ->
        outputPrefix = "${params.outputDir}/allele_minima/${id}"
        inferred_tree.moveTo("${outputPrefix}_inferred_tree.txt")
        inferred_results.moveTo("${outputPrefix}_inferred_results.json")
        timing.moveTo("${outputPrefix}_timing.txt")
    }

    // run Pairtree
    // simulation | create_pairtree_input | pairtree | map { results, best_tree, timing, id ->
        // outputPrefix = "${params.outputDir}/pairtree/${id}"
        // results.moveTo("${outputPrefix}_results.npz")
        // best_tree.moveTo("${outputPrefix}_best_tree.txt")
        // timing.moveTo("${outputPrefix}_timing.txt")
    // }
}
