process AGGREGATE_QUERY {
    label 'process_medium'
    conda '/home/rschwartz/anaconda3/envs/scanpyenv'

    input:
    tuple val(study), path(query_paths)

    output:
    path "${study}_combined.h5ad", emit: combined_query_adata

    script:
    """
    python $projectDir/bin/aggregate_query.py \\
        --query_paths ${query_paths} \\
        --study ${study} \\
        --n_cells ${params.subsample_combined} \\
        --seed ${params.seed}
    """
}
