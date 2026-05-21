# Refactoring `get_census` in `utils.py`

## Problem

`get_census()` in `nextflow_eval_pipeline/bin/utils.py` conflates three distinct operations into one function:

1. **Filter obs metadata** from CellXGene by organism, organ, disease, collections
2. **Split and extract** — subsample and fetch embeddings independently per `split_column` value (e.g., per `dataset_title`)
3. **Whole extraction** — subsample and fetch embeddings for all filtered cells together

This causes the Census API to be called **twice** for the same data: once inside `split_and_extract_data` (N times, once per split), and once again at the bottom of `get_census` to build `refs["whole cortex"]`. For large references this is slow and wasteful.

The calling code in `get_census_adata.py` then saves all split outputs as separate h5ad files, and a separate `combine_refs.py` step re-concatenates them — reconstructing something that was already available as `refs["whole cortex"]`.

---

## Proposed Refactor

Break `get_census` into three focused functions:

### 1. `get_filtered_obs(census, organism, organ, collections, original_celltypes) -> pd.DataFrame`

Already partially exists as `get_filtered_obs()` but without collection filtering or author label mapping. Extend it to:
- Filter by organism, organ, disease, primary data
- Join with `dataset_info` to get collection metadata
- Filter to `ref_collections`
- Optionally map author labels via `map_author_labels()`

Returns a plain obs DataFrame — no data download yet.

### 2. `extract_ref(obs, census, organism, subsample, relabel_path, ref_keys, seed) -> AnnData`

Replaces the current `extract_data()`. Takes a pre-filtered obs DataFrame and:
- Subsamples cell IDs via `subsample_cells()`
- Fetches `get_anndata()` from Census with `obs_embeddings=["scvi"]`
- Runs gene filtering (`filter_genes`)
- Merges dataset_info back onto obs
- Applies `relabel()`

Returns a single AnnData. No split logic, no hidden Census opens.

### 3. `split_refs(obs, split_column, ...) -> dict[str, AnnData]`  *(optional, for per-dataset use cases)*

Thin wrapper that calls `extract_ref()` per unique value in `split_column`. Keep this separate so callers that need per-dataset objects can use it explicitly, and callers that just want the combined reference don't pay for it.

---

## Resulting call patterns

**Combined reference only (this pipeline):**
```python
census = cellxgene_census.open_soma(census_version=census_version)
obs    = get_filtered_obs(census, organism, organ, ref_collections, original_celltypes)
ref    = extract_ref(obs, census, organism, subsample, relabel_path, ref_keys, seed)
ref.write_h5ad("whole_cortex.h5ad")
```

**Per-dataset references (eval pipeline):**
```python
census = cellxgene_census.open_soma(census_version=census_version)
obs    = get_filtered_obs(census, organism, organ, ref_collections, original_celltypes)
refs   = split_refs(obs, split_column="dataset_title", census=census, ...)
for name, adata in refs.items():
    adata.write_h5ad(f"refs/{name}.h5ad")
```

---

## What to remove

| Current | Status |
|---|---|
| `get_census()` | Remove — replaced by direct composition of the three functions above |
| `split_and_extract_data()` | Fold into `split_refs()` |
| `extract_data()` | Rename/simplify to `extract_ref()`, remove `dims` param (unused) |
| `combine_refs.py` + `COMBINE_REFS` module | Remove — no longer needed |
| `GET_CENSUS_ADATA` module | Replace with `GET_COMBINED_REF` (already done in this pipeline) |

## Notes

- `dims` is passed through `get_census` → `split_and_extract_data` → `extract_data` but is never used inside `extract_data` — can be deleted entirely.
- `cell_columns` is hardcoded inside `get_census` but passed as a parameter to `extract_data` — should just be a constant or default in `extract_ref`.
- The Census connection (`census = cellxgene_census.open_soma(...)`) is opened inside `get_census` and never explicitly closed. Refactored code should use it as a context manager: `with cellxgene_census.open_soma(...) as census:`.
