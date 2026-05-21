process PROCESS_REF_SEURAT {
    label 'process_high'
    conda '/home/rschwartz/anaconda3/envs/r4.3'

    input:
    path h5ad_file

    output:
    path "${h5ad_file.getName().replace('.h5ad', '.rds')}", emit: ref_rds

    script:
    """
    Rscript $projectDir/bin/ref_preprocessing.R \\
        --h5ad_file          ${h5ad_file} \\
        --normalization_method ${params.normalization_method} \\
        --dims               ${params.seurat_dims} \\
        --nfeatures          ${params.nfeatures} \\
        ${params.batch_correct ? '--batch_correct' : ''}
    """
}
