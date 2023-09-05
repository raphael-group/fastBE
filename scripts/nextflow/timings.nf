params.proj_dir      = "/n/fs/ragr-research/projects/vafpp"
params.sim_script    = "${params.proj_dir}/scripts/simulation.py"
params.python_script = "${params.proj_dir}/scripts/vafpp_lp.py"
params.cpp_command   = "${params.proj_dir}/build/src/vafpp"

params.nmutations = [500, 1000, 2000, 4000]
params.nclones    = [250, 500, 750, 1000]
params.nsamples   = [50, 100, 200, 500]
params.seeds      = [0, 1, 2, 3, 4, 5]
params.coverage   = [50] 

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

process regress_python {
    cpus 4
    memory '32 GB'
    time '59m'
    errorStrategy 'ignore'

    input:
        tuple path(clonal_matrix), path(mut_clone_mapping), path(freq_matrix), path(total_matrix),
              path(clone_tree), path(usage_matrix), path(variant_matrix), val(id)

    output:
        tuple file("python_results.json"), val(id)

    // module load gurobi
    """
    python '${params.python_script}' --tree ${clone_tree} --frequency-matrix ${freq_matrix} --output python
    """
}

process regress_cpp {
    cpus 1
    memory '2 GB'
    time '59m'

    input:
        tuple path(clonal_matrix), path(mut_clone_mapping), path(freq_matrix), path(total_matrix),
              path(clone_tree), path(usage_matrix), path(variant_matrix), val(id)

    output:
        tuple file("cpp_results.json"), val(id)

    """
    '${params.cpp_command}' regress ${clone_tree} ${freq_matrix} --output cpp
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

    simulation =  parameter_channel | create_sim 
    cpp_res = simulation | regress_cpp
    py_res  = simulation | regress_python

    cpp_res | map { result, name ->
      result.moveTo("nextflow_results/timings/${name}_cpp_results.json")
    }

    py_res | map { result, name ->
      result.moveTo("nextflow_results/timings/${name}_python_results.json")
    }
}
