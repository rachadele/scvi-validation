# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Purpose

Validation of a scVI-based cell type annotation pipeline. Assesses whether scVI embeddings and RF/kNN classifiers produce meaningful, trustworthy predictions.

## Running the Pipeline

```bash
# Human
nextflow run main.nf -params-file params.hs.json -profile conda

# Mouse
nextflow run main.nf -params-file params.mm.json -profile conda

# With per-dataset subsampling (N cells per dataset_key group per query file)
nextflow run main.nf -params-file params.hs.json --subsample_per_dataset 50 --dataset_key sample_id
```

## Upstream Pipeline

This pipeline is a sibling to `../annotation-benchmark`. It reuses that pipeline's Python scripts (via `params.eval_pipeline_dir`) rather than duplicating them. The MAP_QUERY, GET_CENSUS_ADATA, and RUN_SETUP modules all call scripts from `${params.eval_pipeline_dir}/bin/`.

The only new scripts live in `bin/` here:
- `embedding_coherence.py` — centroid cosine similarity (ref vs query) per cell type
- `joint_umap.py` — joint UMAP of ref + query on `obsm["scvi"]`
- `compare_classifiers.py` — RF vs weighted kNN F1 comparison; imports `rfc_pred` from eval pipeline's `utils.py` via PYTHONPATH

## Pipeline Inputs → Outputs

```
params.queries_adata  ─┐
params.relabel_q      ─┤→ MAP_QUERY (process_query_val.py) → *_processed.h5ad
params.relabel_r       │
CellXGene Census      ─┘→ PREPARE_REFERENCES → ref.h5ad

(processed_query, ref) pairs → VALIDATION_PIPELINE:
  ├── EMBEDDING_COHERENCE → results/embedding_coherence/{study}/{ref}/*.coherence.{tsv,png}
  ├── JOINT_UMAP          → results/umap/{study}/{ref}/*.umap.png
  └── COMPARE_CLASSIFIERS → results/classifier_comparison/{study}/{ref}/{query}/*.{tsv,png}
```

## Key Parameters

| param | default | description |
|---|---|---|
| `subsample_query` | 100 | cells sampled per query file |
| `n_neighbors` | 15 | kNN k (also used for UMAP graph) |
| `eval_pipeline_dir` | `/space/grp/rschwartz/rschwartz/annotation-benchmark` | path to sibling pipeline |

## Cell Type Hierarchy

Labels are structured as `subclass → class → family → global` (controlled by `ref_keys`). This hierarchy is defined in the Census maps:
- Human: `meta/census_map_human.tsv` (in eval pipeline)
- Mouse: `meta/census_map_mouse_author.tsv` (in eval pipeline)
