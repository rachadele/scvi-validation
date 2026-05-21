#!/usr/bin/env python3
"""
embedding_coherence.py

Computes pairwise cosine similarity between all ref cell type centroids and all
query cell type centroids (full cross matrix, not just diagonal). High diagonal
values indicate the embedding aligns same cell types across ref and query;
low off-diagonal values indicate separation between types.

Outputs (per ref_key):
  {query_name}_{ref_name}.{key}.coherence.tsv  — matrix (ref types × query types)
  {query_name}_{ref_name}.{key}.coherence.png  — heatmap
"""

import argparse
import os
import sys
import numpy as np
import pandas as pd
import anndata as ad
from sklearn.metrics.pairwise import cosine_similarity


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('--ref_path',   type=str, required=True)
    parser.add_argument('--query_path', type=str, required=True)
    parser.add_argument('--ref_keys',   type=str, nargs='+',
                        default=['subclass', 'class', 'family', 'global'])
    parser.add_argument('--seed', type=int, default=42)
    if __name__ == '__main__':
        known_args, _ = parser.parse_known_args()
        return known_args
    return parser.parse_args()


def compute_centroids(adata, label_col, embedding_key='scvi'):
    """Return (labels, centroid_matrix) sorted by label."""
    X = adata.obsm[embedding_key]
    labels = adata.obs[label_col].astype(str)
    unique_labels = sorted(labels.unique())
    centroids = np.vstack([X[labels == lbl].mean(axis=0) for lbl in unique_labels])
    return unique_labels, centroids


def main():
    args = parse_arguments()
    ref_name   = os.path.basename(args.ref_path).replace('.h5ad', '')
    query_name = os.path.basename(args.query_path).replace('.h5ad', '')
    prefix     = f'{query_name}_{ref_name}'

    ref   = ad.read_h5ad(args.ref_path,   backed='r')
    query = ad.read_h5ad(args.query_path, backed='r')

    for key in args.ref_keys:
        if key not in ref.obs.columns or key not in query.obs.columns:
            print(f'Skipping {key}: not in both ref and query obs', file=sys.stderr)
            continue

        ref_labels,   ref_centroids   = compute_centroids(ref,   key)
        query_labels, query_centroids = compute_centroids(query, key)

        # Full pairwise cosine similarity: rows = ref types, cols = query types
        sim_matrix = cosine_similarity(ref_centroids, query_centroids)

        df = pd.DataFrame(sim_matrix, index=ref_labels, columns=query_labels)
        df.index.name   = 'ref_cell_type'
        df.columns.name = 'query_cell_type'
        df.to_csv(f'{prefix}.{key}.coherence.tsv', sep='\t')


if __name__ == '__main__':
    main()
