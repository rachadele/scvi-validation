process COMPARE_COHERENCE_METHODS {
    label 'process_low'
    conda '/home/rschwartz/anaconda3/envs/scanpyenv'
    publishDir "${params.outdir}/embedding_coherence_plots/comparison", mode: 'copy'

    input:
    path scvi_tsvs
    path seurat_tsvs
    val  ref_keys

    output:
    path "comparison/*.tsv", emit: comparison_tsv
    path "comparison/*.png", emit: comparison_plots

    script:
    """
    python $projectDir/bin/compare_coherence_methods.py \\
        --scvi_tsvs   ${scvi_tsvs}   \\
        --seurat_tsvs ${seurat_tsvs} \\
        --ref_keys    ${ref_keys}
    """
}
