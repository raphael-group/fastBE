params.proj_dir      = "/n/fs/ragr-research/projects/fastBE"
params.outputDir     = "/n/fs/ragr-research/projects/fastBE/nextflow_results/clustering/"
params.sim_script    = "${params.proj_dir}/scripts/simulation.py"

params.pyclonevi = "/n/fs/ragr-data/users/schmidt/miniconda3/envs/pyclone-vi/bin/pyclone-vi"
params.fastbe = "${params.proj_dir}/build/src/fastbe"

params.nmutations = [100]//, 250, 500]
params.nclones    = [5, 10, 20]
params.nsamples   = [3, 5, 10, 20]
params.seeds      = [0, 1, 2, 3, 4, 5]
params.coverage   = [20]

process create_sim {
    cpus 1
    memory '4 GB'
    time '59m'
    errorStrategy 'ignore'

    publishDir "${params.outputDir}/ground_truth/m${mutations}_n${clones}_s${samples}_c${coverage}_r${seed}", mode: 'copy', overwrite: true

    input:
        tuple val(mutations), val(samples), val(clones), val(coverage), val(seed)

    output:
        tuple file("sim_clonal_matrix.txt"), file("sim_mutation_to_clone_mapping.txt"), file("sim_obs_frequency_matrix.txt"), 
              file("sim_obs_full_frequency_matrix.txt"), file("sim_total_matrix.txt"), file("sim_tree.txt"), 
              file("sim_usage_matrix.txt"), file("sim_variant_matrix.txt"),
              val(clones), val("m${mutations}_n${clones}_s${samples}_c${coverage}_r${seed}")

    """
    python '${params.sim_script}' --mutations ${mutations} --samples ${samples} --clones ${clones} --coverage ${coverage} --seed ${seed} --output sim
    """
}

process fastbe {
    cpus 16
    memory '4 GB'
    time '48h'
    errorStrategy 'ignore'

    publishDir "${params.outputDir}/fastbe/${id}", mode: 'copy', overwrite: true

    input:
        tuple path(clonal_matrix), path(mut_clone_mapping), path(collapsed_freq_matrix), path(freq_matrix), path(total_matrix),
              path(clone_tree), path(usage_matrix), path(variant_matrix), val(clones), val(id)

    output:
        tuple file("inferred_tree.txt"), file("inferred_results.json"), 
              file("inferred_clustering.csv"), file("inferred_clustering_results.json"),
              file("inference_timing.txt"), file("clustering_timing.txt"), val(id)

    """
    /usr/bin/time -v '${params.fastbe}' search ${freq_matrix} -o inferred 2>> inference_timing.txt
    /usr/bin/time -v '${params.fastbe}' cluster inferred_tree.txt ${freq_matrix} -k ${clones} -o inferred 2>> clustering_timing.txt
    """
}

process pyclonevi {
    cpus 16
    memory '4 GB'
    time '48h'

    publishDir "${params.outputDir}/pyclonevi/${id}", mode: 'copy', overwrite: true

    input:
        tuple path(clonal_matrix), path(mut_clone_mapping), path(collapsed_freq_matrix), path(freq_matrix), path(total_matrix),
              path(clone_tree), path(usage_matrix), path(variant_matrix), val(clones), val(id)

    output:
        tuple file("pyclonevi_input.tsv"), file("results.tsv"), file("results.csv"), file("timing.txt")

    """
    python '${params.proj_dir}/scripts/processing/make_pyclonevi_input.py' -v ${variant_matrix} -t ${total_matrix} > pyclonevi_input.tsv
    /usr/bin/time -v '${params.pyclonevi}' fit -i pyclonevi_input.tsv -o output.h5 -c ${clones} 2>> timing.txt
    /usr/bin/time -v '${params.pyclonevi}' write-results-file -i output.h5 -o results.tsv
    python '${params.proj_dir}/scripts/processing/parse_pyclonevi_results.py' results.tsv > results.csv
    """
}

workflow {
    parameter_channel = channel.fromList(params.nmutations)
                               .combine(channel.fromList(params.nsamples))
                               .combine(channel.fromList(params.nclones))
                               .combine(channel.fromList(params.coverage))
                               .combine(channel.fromList(params.seeds))

    simulation = parameter_channel | create_sim 
    simulation | pyclonevi 
}
