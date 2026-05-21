"""
coherence_utils.py — shared helpers for embedding coherence scripts.
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns


def plot_heatmap(df, title, out_path, vmin=-1, vmax=1, cmap='coolwarm',
                 linewidths=0.3, scale=0.45, base_size=8):
    """Render a cosine-similarity heatmap and save to out_path."""
    fig_h = max(base_size, len(df.index)   * scale)
    fig_w = max(base_size, len(df.columns) * scale)
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    annot = len(df.index) <= 20 and len(df.columns) <= 20
    sns.heatmap(df, ax=ax, vmin=vmin, vmax=vmax, cmap=cmap,
                xticklabels=True, yticklabels=True,
                linewidths=linewidths, linecolor='lightgrey',
                annot=annot, fmt='.2f', annot_kws={'size': 9})
    ax.set_title(title, fontsize=14)
    ax.set_xlabel('Query cell type', fontsize=13)
    ax.set_ylabel('Ref cell type',   fontsize=13)
    plt.xticks(rotation=90, fontsize=11)
    plt.yticks(rotation=0,  fontsize=11)
    ax.collections[0].colorbar.ax.tick_params(labelsize=11)
    plt.tight_layout()
    plt.savefig(out_path, dpi=150)
    plt.close()


def diagonal_series(df):
    """Return a Series of diagonal (same-type) cosine similarities."""
    shared = sorted(set(df.index) & set(df.columns))
    return pd.Series({ct: df.loc[ct, ct] for ct in shared},
                     name='cosine_similarity')


def mean_matrix(matrices):
    """Element-wise mean across study matrices, aligned on the union of labels."""
    all_ref_types   = sorted(set().union(*[set(m.index)   for m in matrices.values()]))
    all_query_types = sorted(set().union(*[set(m.columns) for m in matrices.values()]))
    stack = np.stack([
        df.reindex(index=all_ref_types, columns=all_query_types).values
        for df in matrices.values()
    ])
    mean = np.nanmean(stack, axis=0)
    return pd.DataFrame(mean, index=all_ref_types, columns=all_query_types)
