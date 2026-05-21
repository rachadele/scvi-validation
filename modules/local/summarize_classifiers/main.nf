process SUMMARIZE_CLASSIFIERS {
    label 'process_low'
    conda '/home/rschwartz/anaconda3/envs/scanpyenv'
    publishDir "${params.outdir}/classifier_comparison/summary", mode: 'copy'

    input:
    path comparison_tsvs
    val  ref_keys

    output:
    path "summary.classifier_comparison.tsv",          emit: summary_tsv
    path "summary.classifier_comparison.macro_f1.png", emit: macro_f1_plot
    path "summary.*.perclass_f1.png",                  emit: perclass_plots

    script:
    """
    python $projectDir/bin/summarize_classifiers.py \\
        --comparison_tsvs ${comparison_tsvs} \\
        --ref_keys ${ref_keys}
    """
}
