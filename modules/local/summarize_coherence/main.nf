process SUMMARIZE_COHERENCE {
    label 'process_low'
    conda '/home/rschwartz/anaconda3/envs/scanpyenv'
    publishDir "${params.outdir}/embedding_coherence_plots", mode: 'copy'

    input:
    path coherence_tsvs
    val  ref_keys
    val  method

    output:
    path "${method}/summary/summary.${method}.*.mean_coherence.tsv", emit: summary_tsv
    path "${method}/summary/summary.${method}.*.mean_coherence.png", emit: mean_heatmaps
    path "${method}/summary/summary.${method}.*.diagonal_strip.png", emit: diagonal_strips

    script:
    """
    python $projectDir/bin/summarize_coherence.py \\
        --coherence_tsvs ${coherence_tsvs} \\
        --ref_keys       ${ref_keys} \\
        --method         ${method}
    """
}
