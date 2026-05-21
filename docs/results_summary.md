# scVI validation — results summary (mouse, 2024-07-01)

Reference: `whole_cortex` (Census). Queries: 7 GSE studies subsampled to ~200–1000 cells. Per cell type the pipeline reports (a) centroid cosine similarity between query and reference cells in `obsm["scvi"]` ("embedding coherence") and (b) per-class F1 from RF and weighted kNN classifiers trained on the reference.

The headline: scVI embeddings are largely faithful — most subclass diagonals sit at 0.78–0.95. The classifier failures cluster on the cell types where the off-diagonal coherence is high, i.e. where two reference cell types occupy overlapping embedding space.

## Cell-type pairs with correlated embeddings → degraded classification

Values below are mean centroid cosines from `embedding_coherence_plots/scvi/summary/summary.scvi.subclass.mean_coherence.tsv` and per-class F1 from `classifier_comparison/summary/summary.classifier_comparison.tsv`.

### 1. Microglia ↔ Macrophage (CNS macrophage family)
- Off-diagonal: query Macrophage vs ref Microglia centroid = **0.91**, query Microglia vs ref Macrophage = **0.54** — query Macrophages look *more like* the Microglia centroid than the Macrophage centroid (0.68).
- Microglia self-coherence is high (0.85), but the centroids of the two classes are nearly indistinguishable.
- Classifier collapse:
  - **Microglia RF F1 = 0.00** in GSE247339.1, GSE247339.2, GSE214244.1, GSE185454, GSE181021.2 (5/6 studies with Microglia present).
  - kNN partially rescues it (0.04–0.67) but is unreliable.
  - Macrophage F1 = 0.12–0.25 across studies.
- At the **family** level (CNS macrophage), F1 jumps to 0.93–1.0 — the two classes are only separable when collapsed.

### 2. OPC ↔ Neural stem cell
- Off-diagonal: NSC vs OPC centroid = **0.57**; OPC vs NSC = 0.17. NSC self-coherence is the lowest in the panel (**0.57**).
- Classifier:
  - NSC F1 = 0.21–0.46 across studies (poor in both RF and kNN).
  - OPC has a striking RF/kNN gap: **RF F1 = 0.00** for OPC in GSE247339.1, GSE214244.1, GSE185454; kNN recovers to 0.29–0.86. RF appears to systematically misroute OPCs (likely into Oligodendrocyte / NSC).
- OPC also pulls toward Astrocyte (cosine 0.49) and Oligodendrocyte (0.47), consistent with its glial-progenitor identity.

### 3. Vascular subclasses: Endothelial ↔ Pericyte ↔ SMC ↔ VLMC
- Off-diagonal cosines within the vascular block:
  - Pericyte vs SMC = **0.91**, Pericyte vs Vascular = 0.83, SMC vs Vascular = 0.80, Endothelial vs Pericyte = 0.70, Endothelial vs SMC = 0.70.
- Self-coherence is high for each individually (Endothelial 0.95, Pericyte 0.93, SMC 0.91), so cells are well-clustered — they just cluster *near each other*.
- Classifier impact is modest at subclass level (Endothelial F1 ≈ 0.91–0.99, Pericyte F1 ≈ 0.67) because the dominant subclass tends to win, but fine-grained vascular distinctions are unstable. At the family level (Vascular) F1 = 0.93–1.00.

### 4. Hippocampal excitatory subtypes: CA1-ProS ↔ CA3 ↔ SUB-ProS ↔ L5 ET ↔ DG
- Off-diagonal cosines:
  - CA1-ProS ↔ L5 ET = **0.75**, CA1-ProS ↔ SUB-ProS = **0.77**, CA3 ↔ CA1-ProS = 0.65, CA3 ↔ SUB-ProS = 0.54, CA2-IG-FC ↔ CA1-ProS = 0.60, DG ↔ L2/3-6 IT = 0.52.
- These are biologically related projection-neuron classes from cortex/hippocampus and they share embedding territory.
- Classifier impact in GSE185454 (the hippocampal study):
  - CA1-ProS F1 = 0.57–0.62, CA3 F1 = 0.59–0.65, DG F1 = 0.85–0.95.
  - At the **class** level (Hippocampal neuron) F1 = 0.97 — the within-region splits are what fail, not the broader call.

### 5. GABAergic interneuron subtypes
- The GABAergic *family* self-coherence drops to **0.70** (vs >0.83 for most other families). At subclass level, PVALB↔SST↔VIP↔LAMP5↔SNCG centroids partially overlap (cross-cosines in the 0.3–0.5 range).
- This study had small interneuron support, so per-subclass F1 is noisy, but it's the family with the weakest within-family separation.

## Where scVI vs Seurat diverge

From `embedding_coherence_plots/comparison/comparison/comparison.{class,family,subclass}.tsv`:

- **OPC**: scVI = 0.83, Seurat = **−0.05** — Seurat embeddings fail entirely on OPC self-coherence; scVI is the only viable embedding for this class.
- **Neural stem cell**: scVI = 0.57, Seurat = 0.35 — scVI better but both weak.
- Seurat is modestly better for Astrocyte, hippocampal subclasses (CA1-ProS, CA3, DG), GABAergic family, and Oligodendrocyte (~0.05–0.15 higher coherence). These are the cell types where Seurat's reference-anchor approach wins on tightness.
- Microglia, Macrophage, Cajal-Retzius, Vascular, Endothelial: within ~0.02 of each other.

## Takeaways

1. **scVI embeddings preserve most subclass identity** (median diagonal ≈ 0.83), so embedding-based classifiers are reasonable in principle.
2. **The classification failures are not random** — they map onto a small number of biologically-related embedding-overlap pairs: Microglia/Macrophage, OPC/NSC, vascular subtypes, hippocampal projection-neuron subtypes, GABAergic interneurons. Collapsing to the family level recovers F1 ≥ 0.9 for almost all of these.
3. **RF is brittle on the overlapping classes** (frequent F1 = 0); weighted kNN is more forgiving and should be preferred when those classes matter.
4. **OPC is a special case**: scVI handles it well, Seurat does not — argues for scVI when OPC predictions are needed.
5. Consumers of these predictions should treat the listed pairs as fundamentally ambiguous at subclass resolution and either (a) report at family/class level, or (b) flag predictions among these pairs as low-confidence.
