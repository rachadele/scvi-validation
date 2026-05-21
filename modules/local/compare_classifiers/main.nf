process COMPARE_CLASSIFIERS {
    label 'process_medium'
    conda '/home/rschwartz/anaconda3/envs/scanpyenv'
    publishDir "${params.outdir}/classifier_comparison/${query_name}/${ref_name}", mode: 'copy'

    input:
    tuple val(query_path), val(ref_path)
    val ref_keys

    output:
    path "*classifier_comparison.tsv", emit: comparison_tsv
    path "*classifier_comparison.png", emit: comparison_plot

    script:
    ref_name   = ref_path.getName().replace('.h5ad', '')
    query_name = query_path.getName().replace('.h5ad', '')
    """
    python $projectDir/bin/compare_classifiers.py \\
        --query_path ${query_path} \\
        --ref_path ${ref_path} \\
        --mapping_file ${params.relabel_r} \\
        --ref_keys ${ref_keys} \\
        --seed ${params.seed} \\
        --n_neighbors ${params.n_neighbors}
    """
}
