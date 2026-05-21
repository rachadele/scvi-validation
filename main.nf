#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { PREPARE_REFERENCES  } from "$projectDir/subworkflows/local/prepare_references/main"
include { PREPROCESS          } from "$projectDir/subworkflows/local/prepare_seurat/main"
include { VALIDATION_PIPELINE } from "$projectDir/subworkflows/local/validation_pipeline/main"

workflow {

    query_paths_adata = Channel.fromPath(params.queries_adata)
    relabel_q_paths   = Channel.fromPath(params.relabel_q)

    ref_collections = params.ref_collections.collect { "\"${it}\"" }.join(' ')

    PREPARE_REFERENCES(
        params.organism,
        params.census_version,
        ref_collections
    )

    PREPROCESS(
        query_paths_adata,
        relabel_q_paths,
        PREPARE_REFERENCES.out.ref_adata,
        PREPARE_REFERENCES.out.model_path,
        params.ref_keys.join(' ')
    )

    VALIDATION_PIPELINE(
        PREPROCESS.out.combined_query_adata,
        PREPARE_REFERENCES.out.ref_adata,
        params.ref_keys.join(' '),
        PREPROCESS.out.query_rds,
        PREPROCESS.out.ref_rds
    )
}

workflow.onComplete {
    println(workflow.success ?
        """
        ================================================================================
        Validation pipeline complete
        --------------------------------------------------------------------------------
        Run as      : ${workflow.commandLine}
        Completed at: ${workflow.complete}
        Duration    : ${workflow.duration}
        Output dir  : ${params.outdir}
        ================================================================================
        """.stripIndent() :
        "Failed: ${workflow.errorReport}\nExit status: ${workflow.exitStatus}"
    )
}

workflow.onError = {
    println "Error: check '.nextflow.log'"
}
