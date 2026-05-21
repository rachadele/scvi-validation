#!/usr/bin/env python3
"""
summarize_coherence.py

Collects per-study coherence TSVs for a single method (scvi or seurat) and produces:
  1. Mean cosine similarity heatmap across all studies (per ref_key)
  2. Strip plot of diagonal (same cell type) similarities per study

Outputs (under {method}/summary/):
  summary.{method}.{key}.mean_coherence.tsv
  summary.{method}.{key}.mean_coherence.png
  summary.{method}.{key}.diagonal_strip.png
"""

import argparse
import os
import sys
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from coherence_utils import plot_heatmap, diagonal_series, mean_matrix


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('--coherence_tsvs', type=str, nargs='+', required=True)
    parser.add_argument('--ref_keys', type=str, nargs='+',
                        default=['subclass', 'class', 'family', 'global'])
    parser.add_argument('--method', type=str, choices=['scvi', 'seurat'], required=True)
    if __name__ == '__main__':
        known_args, _ = parser.parse_known_args()
        return known_args
    return parser.parse_args()


def tsv_suffix(method):
    return '.seurat_coherence.tsv' if method == 'seurat' else '.coherence.tsv'


def load_matrices(tsvs, key, method):
    """Load all TSVs for a given key and method, return dict of {study: DataFrame}."""
    suffix = f'.{key}{tsv_suffix(method)}'
    matrices = {}
    for path in tsvs:
        if not os.path.basename(path).endswith(suffix):
            continue
        df = pd.read_csv(path, sep='\t', index_col=0)
        df.index.name = 'ref_cell_type'
        study = os.path.basename(path).split('_combined')[0]
        matrices[study] = df
    return matrices


def plot_diagonal_strip(matrices, key, method, outdir):
    """Strip plot: one point per (study, cell_type) for diagonal (same-type) similarity."""
    records = []
    for study, df in matrices.items():
        for ct, val in diagonal_series(df).items():
            records.append({'study': study, 'cell_type': ct, 'cosine_similarity': val})
    if not records:
        return

    strip_df = pd.DataFrame(records)
    order = (strip_df.groupby('cell_type')['cosine_similarity']
             .median().sort_values(ascending=False).index.tolist())

    fig_w = max(10, len(order) * 0.55)
    fig, ax = plt.subplots(figsize=(fig_w, 8))
    sns.stripplot(data=strip_df, x='cell_type', y='cosine_similarity', order=order,
                  hue='study', ax=ax, jitter=True, size=7, alpha=0.8)
    ax.axhline(0, color='grey', linewidth=0.8, linestyle='--')
    ax.set_title(f'Diagonal cosine similarity by cell type [{method.upper()}] ({key})',
                 fontsize=14)
    ax.set_xlabel('Cell type', fontsize=13)
    ax.set_ylabel('Cosine similarity (same-type ref vs query centroid)', fontsize=13)
    plt.xticks(rotation=90, fontsize=11)
    plt.yticks(fontsize=11)
    ax.legend(title='Study', bbox_to_anchor=(1.01, 1), loc='upper left',
              fontsize=10, title_fontsize=11)
    plt.tight_layout()
    plt.savefig(os.path.join(outdir, f'summary.{method}.{key}.diagonal_strip.png'), dpi=150)
    plt.close()


def main():
    args = parse_arguments()

    outdir = os.path.join(args.method, 'summary')
    os.makedirs(outdir, exist_ok=True)

    for key in args.ref_keys:
        matrices = load_matrices(args.coherence_tsvs, key, args.method)
        if not matrices:
            print(f'No TSVs found for key={key} method={args.method}, skipping')
            continue

        mean_df = mean_matrix(matrices)
        mean_df.to_csv(
            os.path.join(outdir, f'summary.{args.method}.{key}.mean_coherence.tsv'), sep='\t')

        title = (f'Mean centroid cosine similarity [{args.method.upper()}] ({key})\n'
                 f'n = {len(matrices)} studies')
        plot_heatmap(mean_df, title,
                     os.path.join(outdir, f'summary.{args.method}.{key}.mean_coherence.png'),
                     scale=0.55, base_size=10)

        plot_diagonal_strip(matrices, key, args.method, outdir)


if __name__ == '__main__':
    main()
