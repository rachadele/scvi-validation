process COMBINE_REFS {
    label 'process_medium'
    conda '/home/rschwartz/anaconda3/envs/scanpyenv'

    input:
    path ref_paths

    output:
    path "whole_cortex.h5ad", emit: combined_ref

    script:
    """
    python $projectDir/bin/combine_refs.py \\
        --ref_paths ${ref_paths} \\
        --output whole_cortex.h5ad
    """
}
