nextflow.enable.dsl=2

params.root_dir = '/n/fs/ragr-research/projects/vafpp/'
params.ground_truth_dir = 'nextflow_results/ground_truth/'
params.algorithms = [
    ['pairtree', 'nextflow_results/pairtree/', '_best_tree.txt'],
    ['allele_minima', 'nextflow_results/allele_minima/', '_inferred_tree.txt']
]

def trimSuffix(original, suffix) {
	if(original.endsWith(suffix)) {
	    return original.substring(0, original.length() - suffix.length())
	}
	return original
}

process EvaluateTrees {
    errorStrategy 'ignore'

    input:
        tuple val(algo), val(name), path(ground_truth), path(inferred_tree)

    output:
        tuple path("result.json"), val(algo), val(name)

    """
    python ${params.root_dir}/scripts/score_result.py $ground_truth $inferred_tree > result.json
    """
}

workflow {
    ground_truth_trees_ch = Channel
        .fromPath(params.ground_truth_dir + '*_tree.txt')
        .map { file -> trimSuffix(file.baseName, '_tree') }

    algorithms_ch = Channel.fromList(params.algorithms)
    eval_channel = algorithms_ch.combine(ground_truth_trees_ch).map { 
        algo = it[0]
        ground_truth_tree = "${params.root_dir}/${params.ground_truth_dir}/${it[3]}_tree.txt"
        inferred_tree = "${params.root_dir}/${it[1]}/${it[3]}${it[2]}"
        [algo, it[3], ground_truth_tree, inferred_tree]
    }

    eval_channel | EvaluateTrees | map { result, algo, name ->
        result.moveTo("nextflow_results/evaluation/${algo}_${name}.json")
    }
}

