params.proj_dir = "/n/fs/ragr-research/projects/vafpp"
params.scripts_dir = "${params.proj_dir}/realdata/scripts"
params.output_dir = "${params.proj_dir}/realdata"

params.pairtree_bin  = "${params.proj_dir}/dependencies/pairtree/bin/pairtree"
params.pairtree_parse_output = "${params.proj_dir}/scripts/parse_pairtree_output.py"                                                                                                   

params.allele_minima = "${params.proj_dir}/build/src/vafpp"                                                                                                                            

process parse_pairtree {
  input:
    tuple path(params_json), path(ssm), val(id)
  output:
    tuple file("exp_mutation_clone_mapping.csv"), file("exp_frequency_matrix.txt"), file("exp_variant_matrix.txt"), 
          file("exp_total_matrix.txt"), file("exp_samples.txt"), file("exp_mutations.txt"), val(id)

  """
  python '${params.scripts_dir}/pairtree_input_parse.py' ${ssm} ${params_json} --output exp
  """
}

process allele_minima {                                                                                                                                                                
    cpus 16                                                                                                                                                                            
    memory '2 GB'                                                                                                                                                                      
    time '24h'                                                                                                                                                                         
                                                                                                                                                                                       
    input:                                                                                                                                                                             
        tuple path(mut_clone_mapping), path(freq_matrix), path(variant_matrix),                                                                                     
              path(total_matrix), val(samples), val(mutations), val(id)                                                                                         
                                                                                                                                                                                       
    output:                                                                                                                                                                            
        tuple file("inferred_tree.txt"), file("inferred_results.json"), file("timing.txt"), val(id)                                                                                    
                                                                                                                                                                                       
    """                                                                                                                                                                                
    /usr/bin/time -v '${params.allele_minima}' search ${freq_matrix} -a 1 -s 5000 --output inferred -t ${task.cpus} 2>> timing.txt                                                      
    """                                                                                                                                                                                
}

process pairtree {                                                                                                                                                            [80/1902]
    cpus 16               
    memory '8 GB'                                                                          
    time '24h'                                                                             
                                                                                           
    input:                                                                                 
        tuple file(params_json), file(ssm), val(id)                                        
                                                                                           
    output:                                                                                
        tuple file("results.npz"), file("best_tree.txt"), file("timing.txt"), val(id)      
                                                                                                                                                                                       
    """                                                                                                                                                                                
    /usr/bin/time -v '${params.pairtree_bin}' --parallel 16 --params ${params_json} ${ssm} results.npz 2>> timing.txt
    python '${params.pairtree_parse_output}' results.npz --output best_tree.txt 
    """                                                                                    
}

workflow {
  inputChannel = Channel.fromFilePairs("${params.proj_dir}/realdata/input/*.{ssm,params.json}") 
  // inputChannel | map {[it[1][0], it[1][1], it[0]]} | map {print it}
  inputChannel | map {[it[1][0], it[1][1], it[0]]} | parse_pairtree | allele_minima | map { inferred_tree, inferred_results, timing, id ->                                                                                                  
        output_prefix = "${params.output_dir}/allele_minima/${id}"                                                                                                                       
        inferred_tree.moveTo("${output_prefix}_inferred_tree.txt")                                                                                                                      
        inferred_results.moveTo("${output_prefix}_inferred_results.json")                                                                                                               
        timing.moveTo("${output_prefix}_timing.txt")  
  }

  inputChannel | map {[it[1][0], it[1][1], it[0]]} | pairtree | map { inferred_results, inferred_tree, timing, id ->                                                                                                  
        output_prefix = "${params.output_dir}/pairtree/${id}"                                                                                                                       
        inferred_tree.moveTo("${output_prefix}_best_tree.txt")                                                                                                                      
        inferred_results.moveTo("${output_prefix}_results.npz")                                                                                                               
        timing.moveTo("${output_prefix}_timing.txt")  
  }

}
