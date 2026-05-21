process GET_COMBINED_REF {
    label 'process_high'
    conda '/home/rschwartz/anaconda3/envs/scanpyenv'
    publishDir "${params.outdir}/refs", mode: 'copy'

    input:
    val ref_collections

    output:
    path "whole_cortex.h5ad", emit: combined_ref

    script:
    ref_keys = params.ref_keys.join(' ')
    """
    python $projectDir/bin/get_combined_ref.py \\
        --organism        ${params.organism} \\
        --census_version  ${params.census_version} \\
        --subsample_ref   ${params.subsample_ref} \\
        --relabel_path    ${params.relabel_r} \\
        --ref_collections ${ref_collections} \\
        --ref_keys        ${ref_keys} \\
        --organ           ${params.organ} \\
        --seed            ${params.seed}
    """
}
