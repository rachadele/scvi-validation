process PROCESS_QUERY_SEURAT {
    label 'process_medium'
    conda '/home/rschwartz/anaconda3/envs/r4.3'

    input:
    path h5ad_file

    output:
    path "${h5ad_file.getName().replace('.h5ad', '.rds')}", emit: query_rds

    script:
    """
    Rscript $projectDir/bin/seurat_preprocessing.R \\
        --h5ad_file          ${h5ad_file} \\
        --normalization_method ${params.normalization_method} \\
        --dims               ${params.seurat_dims} \\
        --nfeatures          ${params.nfeatures}
    """
}
