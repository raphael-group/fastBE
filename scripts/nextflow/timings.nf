params.sim_script = '/home/schmidt73/Desktop/VAFPP Optimization/scripts/simulation.py'
params.python_script = '/home/schmidt73/Desktop/VAFPP Optimization/scripts/vafpp_lp.py'
params.cpp_command = '/home/schmidt73/Desktop/VAFPP Optimization/build/src/vafpp'

process create_sim {
    input:
        tuple val(mutations), val(samples), val(clones), val(coverage), val(seed)

    output:
        tuple file("sim_clonal_matrix.txt"), file("sim_mutation_to_clone_mapping.txt"), file("sim_obs_frequency_matrix.txt"), 
              file("sim_total_matrix.txt"), file("sim_tree.txt"), file("sim_usage_matrix.txt"), file("sim_variant_matrix.txt")

    """
    python '${params.sim_script}' --mutations ${mutations} --samples ${samples} --clones ${clones} --coverage ${coverage} --seed ${seed} --output sim
    """
}

process regress_python {
    input:
        tuple path(clonal_matrix), path(mut_clone_mapping), path(freq_matrix), path(total_matrix),
              path(clone_tree), path(usage_matrix), path(variant_matrix)

    output:
        file("python_results.json")

    """
    python '${params.python_script}' --tree ${clone_tree} --frequency-matrix ${freq_matrix} --output python
    """
}

process regress_cpp {
    input:
        tuple path(clonal_matrix), path(mut_clone_mapping), path(freq_matrix), path(total_matrix),
              path(clone_tree), path(usage_matrix), path(variant_matrix)

    output:
        file("cpp_results.json")

    """
    '${params.cpp_command}' regress ${clone_tree} ${freq_matrix} --output cpp
    """
}

workflow {
    simulation = channel.of([50, 10, 25, 50, 0]) | create_sim 
    cpp_res = simulation | regress_cpp
    py_res  = simulation | regress_python
}
