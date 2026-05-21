include { RUN_SETUP       } from "$projectDir/modules/local/run_setup/main"
include { GET_COMBINED_REF } from "$projectDir/modules/local/get_combined_ref/main"

workflow PREPARE_REFERENCES {
    take:
    organism
    census_version
    ref_collections

    main:
    model_path = RUN_SETUP(organism, census_version)
    GET_COMBINED_REF(ref_collections)

    emit:
    model_path = model_path
    ref_adata  = GET_COMBINED_REF.out.combined_ref
}
