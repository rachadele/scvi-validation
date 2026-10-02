# Track 1: Validation of Current Framework

## Goal
Assess whether the current scVI embeddings and RF classifier are producing meaningful, trustworthy predictions before investing in methodological changes.

## Pipeline Context

The current pipeline (now `annotation-benchmark`) works as follows:

1. **Model download** (`bin/setup.py` + `utils.setup()`): Downloads a pre-trained scVI model from CellXGene Census for the given organism and census version.
2. **Reference preparation** (`bin/get_census_adata.py`): Fetches reference cells from Census, stores them as h5ad. Reference embeddings are pre-computed and stored in `adata.obsm["scvi"]`.
3. **Query mapping** (`bin/process_query.py` + `utils.map_query()`): Maps query cells through the pre-trained model using `scvi.model.SCVI.load_query_data()`, then stores latent embeddings in `query.obsm["scvi"]`.
4. **RF prediction** (`bin/predict_scvi.py` + `utils.rfc_pred()`): Trains a `RandomForestClassifier(class_weight='balanced')` on `ref.obsm["scvi"]` and predicts on `query.obsm["scvi"]`. Predictions are aggregated up the cell type hierarchy (subclass → class → family → global).
5. **Metrics** (`bin/classify_all.py`): Computes per-cell-type F1 scores and confusion matrices.

### Input data paths

**Human queries**: `/space/grp/rschwartz/rschwartz/get_gemma_data.nf/all_homo_sapiens_samples/h5ad/*/*h5ad`
**Mouse queries**: `/space/grp/rschwartz/rschwartz/get_gemma_data.nf/study_names_mouse.txt_author_true_process_samples_true/h5ad/**h5ad`

**Relabeling files (query)** — harmonize author cell type labels to a common vocabulary:
- Human: `meta/relabel_homo_sapiens/*_relabel.tsv` (one TSV per study, e.g. `CMC_relabel.tsv`, `Velmeshev-2019.1_relabel.tsv`)
- Mouse: `meta/relabel_mus_musculus/*_relabel.tsv`

**Census maps (reference)** — map Census cell type labels to the hierarchy levels (subclass → class → family → global):
- Human: `meta/census_map_human.tsv` (columns: `cell_type`, `subclass`, `class`, `family`, `global`)
- Mouse: `meta/census_map_mouse_author.tsv` (columns: `author_cell_type`, `subclass`, `class`, `family`, `global`)

**Reference collections used**:
- Human: Allen Brain Cell Atlas neocortex, SEA-AD, primate dlPFC; census version `2024-07-01`
- Mouse: Allen isocortex/hippocampus, motor cortex atlas, Tabula Muris Senis; census version `2024-07-01`

### Key artifacts produced by the pipeline
- `query_mapped.h5ad` — query anndata with `obsm["scvi"]` latent embeddings
- `ref.h5ad` — reference anndata with `obsm["scvi"]` latent embeddings
- TSV files of F1 scores per cell type level
- Confusion matrices

### Key shared code to reuse
- `bin/utils.py`: `map_query()`, `rfc_pred()`, cell type hierarchy aggregation logic
- `bin/classify_all.py`: F1 scoring and confusion matrix logic

## Steps

### 1. Embedding Coherence
- Load paired ref + query h5ads produced by the pipeline
- For each cell type, compute pairwise cosine similarity (or correlation) between reference and query `obsm["scvi"]` embeddings
- Expect: same cell type → high similarity across ref/query; different cell types → low similarity
- Output: heatmap or summary table of within- vs. between-cell-type embedding distances
- **Reuse**: `ref.obsm["scvi"]`, `query.obsm["scvi"]`, cell type labels from `obs`

### 2. Joint UMAP (Reference + Query)
- Concatenate ref and query anndatas (using paired `obsm["scvi"]` embeddings from pipeline output)
- Compute UMAP on the combined latent space
- Color by: cell type label, dataset origin (ref vs. query), prediction correctness
- Output: UMAP plots saved per study/query
- **Reuse**: pipeline-produced h5ads; can use `scanpy.tl.umap` on the `"scvi"` embedding

### 3. Classifier Comparison: RF vs Weighted kNN
- **Contingent on Steps 1 & 2**: if embeddings show good within-cell-type coherence and clear separation in UMAP, kNN is a natural fit (it directly exploits embedding geometry). If embeddings are noisy or poorly separated, kNN may perform worse than RF and the comparison is less meaningful.
- Reuse `utils.rfc_pred()` as the RF baseline
- Implement a weighted kNN classifier on `obsm["scvi"]` using the same ref/query splits
- Use same evaluation metrics as `classify_all.py` (F1 per cell type, macro F1)
- Compare: does kNN outperform RF? Are errors on the same cell types?
- Output: side-by-side F1 comparison table
- **If kNN wins**: replace `rfc_pred()` in `utils.py` with a kNN implementation — this change propagates automatically to the full pipeline and to Track 2 evaluation
- **Reuse**: `utils.rfc_pred()` for baseline; `classify_all.py` logic for scoring

## Deliverables
- Figures for embedding coherence and joint UMAPs
- F1 comparison table (RF vs kNN)
- Short write-up summarizing findings
