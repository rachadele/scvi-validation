#!/usr/bin/env python3
"""
compare_coherence_methods.py

Directly compares mean centroid cosine similarity between scVI and Seurat
across all cell types, per taxonomy level.

Inputs:
  --scvi_tsvs   : summary.{method}.{key}.mean_coherence.tsv from SUMMARIZE_SCVI_COHERENCE
  --seurat_tsvs : summary.{method}.{key}.mean_coherence.tsv from SUMMARIZE_SEURAT_COHERENCE

Outputs per ref_key under comparison/:
  comparison.{key}.selectivity.png   -- diagonal vs max-off-diagonal per cell type, by method
  comparison.{key}.confusion.png     -- confusion profiles for top-N most confusable cell types
  comparison.{key}.tsv               -- long-format table with both scores and difference
"""

import argparse
import os
import re
import sys
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from coherence_utils import diagonal_series


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('--scvi_tsvs',   type=str, nargs='+', required=True)
    parser.add_argument('--seurat_tsvs', type=str, nargs='+', required=True)
    parser.add_argument('--ref_keys',    type=str, nargs='+',
                        default=['subclass', 'class', 'family', 'global'])
    parser.add_argument('--n_confusable', type=int, default=10,
                        help='Number of most confusable cell types for confusion profile plot')
    if __name__ == '__main__':
        known_args, _ = parser.parse_known_args()
        return known_args
    return parser.parse_args()


def find_tsv(tsvs, key):
    for p in tsvs:
        if re.search(rf'summary\.\w+\.{key}\.mean_coherence\.tsv$', os.path.basename(p)):
            return p
    return None


def max_off_diagonal(df):
    """Return a Series of max off-diagonal cosine similarities per row (ref cell type)."""
    shared = sorted(set(df.index) & set(df.columns))
    records = {}
    for ct in shared:
        other_cols = [c for c in df.columns if c != ct]
        records[ct] = df.loc[ct, other_cols].max() if other_cols else float('nan')
    return pd.Series(records, name='max_off_diagonal')



def plot_selectivity_scatter(scvi_df, seurat_df, key, out_path):
    """
    Selectivity scatter: x = diagonal similarity, y = max off-diagonal similarity.
    One point per cell type per method. Points in the lower-right quadrant
    (high diagonal, low off-diagonal) represent cleanly separable cell types.
    Points in the upper-right quadrant are recognisable but promiscuous.
    Points in the upper-left are systematically misclassified.
    """
    scvi_diag   = diagonal_series(scvi_df)
    seurat_diag = diagonal_series(seurat_df)
    scvi_off    = max_off_diagonal(scvi_df)
    seurat_off  = max_off_diagonal(seurat_df)

    fig, ax = plt.subplots(figsize=(9, 8))
    lim = (-0.1, 1.05)

    for diag, off_diag, color, label in [
        (scvi_diag,   scvi_off,   'steelblue', 'scVI'),
        (seurat_diag, seurat_off, 'tomato',    'Seurat'),
    ]:
        common = sorted(set(diag.index) & set(off_diag.index))
        ax.scatter(diag[common], off_diag[common],
                   s=70, alpha=0.8, color=color, edgecolors='white',
                   linewidths=0.5, label=label, zorder=3)
        for ct in common:
            ax.annotate(ct, (diag[ct], off_diag[ct]), fontsize=6, alpha=0.7,
                        xytext=(3, 2), textcoords='offset points')

    # y = x reference: above this line means off-diagonal > diagonal (confused)
    ax.plot(lim, lim, '--', color='grey', linewidth=0.8, alpha=0.5, label='y = x')
    ax.set_xlim(lim)
    ax.set_ylim(lim)
    ax.set_xlabel('Diagonal cosine similarity (same-type)', fontsize=13)
    ax.set_ylabel('Max off-diagonal cosine similarity', fontsize=13)
    ax.set_title(
        f'Selectivity scatter ({key})\n'
        f'lower-right = separable  |  upper-right = promiscuous  |  upper-left = misclassified',
        fontsize=12)
    ax.legend(fontsize=11)
    plt.tight_layout()
    plt.savefig(out_path, dpi=150)
    plt.close()


def _confusable_types(df, n):
    """Top-N cell types by max_off_diagonal / diagonal ratio."""
    diag = diagonal_series(df)
    off  = max_off_diagonal(df)
    shared = sorted(set(diag.index) & set(off.index))
    ratio = off[shared] / diag[shared].replace(0, float('nan'))
    return ratio.nlargest(n).index.tolist()


def plot_confusion_profiles(scvi_df, seurat_df, key, out_path, n_top=10):
    """
    For the top-N most confusable cell types (union across both methods, ranked by
    max-off-diagonal / diagonal ratio), show their full similarity row from each
    method's mean coherence matrix in a side-by-side heatmap.
    """
    confusable = sorted(
        set(_confusable_types(scvi_df,   n_top)) |
        set(_confusable_types(seurat_df, n_top))
    )
    if not confusable:
        return

    shared_query = sorted(set(scvi_df.columns) & set(seurat_df.columns))
    if not shared_query:
        return

    scvi_rows   = scvi_df.reindex(index=confusable, columns=shared_query)
    seurat_rows = seurat_df.reindex(index=confusable, columns=shared_query)

    n_rows = len(confusable)
    n_cols = len(shared_query)
    fig_h  = max(6, n_rows * 0.5 + 2)
    fig_w  = max(16, n_cols * 0.5 + 4)
    annot  = n_rows <= 15 and n_cols <= 15

    fig, axes = plt.subplots(1, 2, figsize=(fig_w, fig_h))
    for ax, data, method_name in zip(axes,
                                     [scvi_rows,   seurat_rows],
                                     ['scVI',       'Seurat']):
        sns.heatmap(data, ax=ax, vmin=-1, vmax=1, cmap='coolwarm',
                    xticklabels=True, yticklabels=True,
                    linewidths=0.3, linecolor='lightgrey',
                    annot=annot, fmt='.2f', annot_kws={'size': 8})
        ax.set_title(f'{method_name}', fontsize=13)
        ax.set_xlabel('Query cell type', fontsize=11)
        ax.set_ylabel('Ref cell type', fontsize=11)
        plt.sca(ax)
        plt.xticks(rotation=90, fontsize=9)
        plt.yticks(rotation=0,  fontsize=9)

    fig.suptitle(
        f'Confusion profiles — top {n_top} most confusable cell types ({key})\n'
        f'ranked by max off-diagonal / diagonal similarity ratio',
        fontsize=13)
    plt.tight_layout()
    plt.savefig(out_path, dpi=150)
    plt.close()


def main():
    args = parse_arguments()
    os.makedirs('comparison', exist_ok=True)

    for key in args.ref_keys:
        scvi_path   = find_tsv(args.scvi_tsvs,   key)
        seurat_path = find_tsv(args.seurat_tsvs, key)

        if scvi_path is None:
            print(f'No scVI summary TSV for key={key}, skipping', file=sys.stderr)
            continue
        if seurat_path is None:
            print(f'No Seurat summary TSV for key={key}, skipping', file=sys.stderr)
            continue

        scvi_df   = pd.read_csv(scvi_path,   sep='\t', index_col=0)
        seurat_df = pd.read_csv(seurat_path, sep='\t', index_col=0)

        scvi_diag   = diagonal_series(scvi_df)
        seurat_diag = diagonal_series(seurat_df)

        common  = sorted(set(scvi_diag.index) & set(seurat_diag.index))
        out_tsv = pd.DataFrame({
            'cell_type':              common,
            'scvi':                   scvi_diag[common].values,
            'seurat':                 seurat_diag[common].values,
            'diff_scvi_minus_seurat': (scvi_diag[common] - seurat_diag[common]).values,
        })
        out_tsv.to_csv(f'comparison/comparison.{key}.tsv', sep='\t', index=False)

        plot_selectivity_scatter(scvi_df, seurat_df, key,
                                 f'comparison/comparison.{key}.selectivity.png')
        plot_confusion_profiles(scvi_df, seurat_df, key,
                                f'comparison/comparison.{key}.confusion.png',
                                n_top=args.n_confusable)
        print(f'Wrote comparison outputs for key={key}')


if __name__ == '__main__':
    main()
