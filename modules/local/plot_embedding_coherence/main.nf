process PLOT_EMBEDDING_COHERENCE {
    label 'process_medium'
    conda '/home/rschwartz/anaconda3/envs/scanpyenv'
    publishDir "${params.outdir}/embedding_coherence_plots", mode: 'copy'

    input:
    path coherence_tsvs

    output:
    path "**/*.png", emit: coherence_plots

    script:
    """
    python $projectDir/bin/plot_embedding_coherence.py \\
        --coherence_tsvs ${coherence_tsvs}
    """
}
