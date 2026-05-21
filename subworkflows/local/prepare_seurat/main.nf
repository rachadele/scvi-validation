include { MAP_QUERY            } from "$projectDir/modules/local/map_query/main"
include { AGGREGATE_QUERY      } from "$projectDir/modules/local/aggregate_query/main"
include { PROCESS_REF_SEURAT   } from "$projectDir/modules/local/process_ref_seurat/main"
include { PROCESS_QUERY_SEURAT } from "$projectDir/modules/local/process_query_seurat/main"

workflow PREPROCESS {
    take:
    query_paths_adata  // path channel: raw query h5ad files
    relabel_q_paths    // path channel: [key, relabel_tsv] tuples
    ref_adata          // path: single whole_cortex.h5ad
    model_path         // val: scVI model path
    ref_keys           // val: space-separated hierarchy keys

    main:
    // Build per-query input tuples for MAP_QUERY
    relabel_q_keyed = relabel_q_paths.map { path ->
        def key = path.getName().split('_relabel.tsv')[0]
        [key, path]
    }

    combined_query_paths = query_paths_adata
        .map { path ->
            def query_name = path.getName().split('.h5ad')[0]
            def query_key  = query_name.split('_')[0]
            [query_key, query_name, path]
        }
        .combine(relabel_q_keyed, by: 0)
        .map { query_key, query_name, query_path, relabel_path ->
            def batch_key = params.batch_keys.get(query_key, "sample_id")
            [query_name, relabel_path, query_path, batch_key]
        }

    MAP_QUERY(model_path, combined_query_paths, ref_keys)

    // Group processed samples by study, then concatenate and subsample
    study_grouped = MAP_QUERY.out.processed_query_adata
        .map { path ->
            def study = path.getName().split('_')[0]
            [study, path]
        }
        .groupTuple()

    AGGREGATE_QUERY(study_grouped)

    // Convert to Seurat RDS for Seurat-based validation
    PROCESS_REF_SEURAT(ref_adata)
    PROCESS_QUERY_SEURAT(AGGREGATE_QUERY.out.combined_query_adata)

    emit:
    combined_query_adata = AGGREGATE_QUERY.out.combined_query_adata
    ref_rds              = PROCESS_REF_SEURAT.out.ref_rds
    query_rds            = PROCESS_QUERY_SEURAT.out.query_rds
}
