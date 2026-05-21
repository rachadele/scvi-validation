process SEURAT_EMBEDDING_COHERENCE {
    label 'process_medium'
    conda '/home/rschwartz/anaconda3/envs/r4.3'
    publishDir "${params.outdir}/seurat_embedding_coherence/${query_name}/${ref_name}", mode: 'copy'

    input:
    tuple val(query_rds), val(ref_rds)
    val ref_keys

    output:
    path "*.seurat_coherence.tsv", emit: coherence_tsv

    script:
    ref_name   = ref_rds.getName().replace('.rds', '')
    query_name = query_rds.getName().replace('.rds', '')
    """
    Rscript $projectDir/bin/seurat_embedding_coherence.R \\
        --query_path ${query_rds} \\
        --ref_path   ${ref_rds}   \\
        --ref_keys   ${ref_keys}  \\
        --dims       ${params.seurat_dims} \\
        --seed       ${params.seed}
    """
}
