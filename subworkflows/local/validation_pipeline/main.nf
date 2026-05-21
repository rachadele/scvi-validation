include { SCVI_EMBEDDING_COHERENCE   } from "$projectDir/modules/local/embedding_coherence/main"
include { SEURAT_EMBEDDING_COHERENCE } from "$projectDir/modules/local/seurat_embedding_coherence/main"
include { PLOT_EMBEDDING_COHERENCE   } from "$projectDir/modules/local/plot_embedding_coherence/main"
include { JOINT_UMAP                 } from "$projectDir/modules/local/joint_umap/main"
include { COMPARE_CLASSIFIERS        } from "$projectDir/modules/local/compare_classifiers/main"
include { SUMMARIZE_COHERENCE as SUMMARIZE_SCVI_COHERENCE   } from "$projectDir/modules/local/summarize_coherence/main"
include { SUMMARIZE_COHERENCE as SUMMARIZE_SEURAT_COHERENCE } from "$projectDir/modules/local/summarize_coherence/main"
include { COMPARE_COHERENCE_METHODS                         } from "$projectDir/modules/local/compare_coherence_methods/main"
include { SUMMARIZE_CLASSIFIERS      } from "$projectDir/modules/local/summarize_classifiers/main"

workflow VALIDATION_PIPELINE {
    take:
    combined_queries   // path channel: one {study}_combined.h5ad per study
    ref_adata          // path: single whole_cortex.h5ad
    ref_keys           // val: space-separated string of hierarchy keys
    query_rds          // path channel: one {study}_combined.rds per study
    ref_rds            // path: single whole_cortex.rds

    main:
    // Pair each study with the single combined reference
    combos = combined_queries.combine(ref_adata)

    SCVI_EMBEDDING_COHERENCE(combos, ref_keys)
    JOINT_UMAP(combos, ref_keys)
    COMPARE_CLASSIFIERS(combos, ref_keys)

    // Seurat embedding coherence
    seurat_combos = query_rds.combine(ref_rds)
    SEURAT_EMBEDDING_COHERENCE(seurat_combos, ref_keys)

    // Plot heatmaps for both scVI and Seurat coherence TSVs
    PLOT_EMBEDDING_COHERENCE(
        SCVI_EMBEDDING_COHERENCE.out.coherence_tsv
            .mix(SEURAT_EMBEDDING_COHERENCE.out.coherence_tsv)
            .collect()
    )

    // Summary across studies — run separately per method
    SUMMARIZE_SCVI_COHERENCE(
        SCVI_EMBEDDING_COHERENCE.out.coherence_tsv.collect(),
        ref_keys,
        'scvi'
    )

    SUMMARIZE_SEURAT_COHERENCE(
        SEURAT_EMBEDDING_COHERENCE.out.coherence_tsv.collect(),
        ref_keys,
        'seurat'
    )

    COMPARE_COHERENCE_METHODS(
        SUMMARIZE_SCVI_COHERENCE.out.summary_tsv.collect(),
        SUMMARIZE_SEURAT_COHERENCE.out.summary_tsv.collect(),
        ref_keys
    )

    SUMMARIZE_CLASSIFIERS(
        COMPARE_CLASSIFIERS.out.comparison_tsv.collect(),
        ref_keys
    )

    emit:
    coherence_tsv           = SCVI_EMBEDDING_COHERENCE.out.coherence_tsv
    seurat_coherence_tsv    = SEURAT_EMBEDDING_COHERENCE.out.coherence_tsv
    coherence_summary       = SUMMARIZE_SCVI_COHERENCE.out.summary_tsv
    umap_plots              = JOINT_UMAP.out.umap_plots
    comparison_tsv          = COMPARE_CLASSIFIERS.out.comparison_tsv
    classifier_summary      = SUMMARIZE_CLASSIFIERS.out.summary_tsv
}
