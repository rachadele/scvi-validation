process SCVI_EMBEDDING_COHERENCE {
    label 'process_medium'
    conda '/home/rschwartz/anaconda3/envs/scanpyenv'
    publishDir "${params.outdir}/embedding_coherence/${query_name}/${ref_name}", mode: 'copy'

    input:
    tuple val(query_path), val(ref_path)
    val ref_keys

    output:
    path "*.coherence.tsv", emit: coherence_tsv

    script:
    ref_name   = ref_path.getName().replace('.h5ad', '')
    query_name = query_path.getName().replace('.h5ad', '')
    """
    python $projectDir/bin/embedding_coherence.py \\
        --query_path ${query_path} \\
        --ref_path ${ref_path} \\
        --ref_keys ${ref_keys} \\
        --seed ${params.seed}
    """
}
