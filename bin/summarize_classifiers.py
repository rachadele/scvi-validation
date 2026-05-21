#!/usr/bin/env python3
"""
summarize_classifiers.py

Collects per-study classifier_comparison.tsv files and produces:
  1. Bar chart of mean macro F1 (RF vs kNN) per hierarchy level
  2. Per-class F1 heatmap (cell type × classifier) averaged across studies, per key

Outputs:
  summary.classifier_comparison.macro_f1.png
  summary.{key}.perclass_f1.png
  summary.classifier_comparison.tsv
"""

import argparse
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('--comparison_tsvs', type=str, nargs='+', required=True)
    parser.add_argument('--ref_keys', type=str, nargs='+',
                        default=['subclass', 'class', 'family', 'global'])
    if __name__ == '__main__':
        known_args, _ = parser.parse_known_args()
        return known_args
    return parser.parse_args()


def main():
    args = parse_arguments()

    dfs = []
    for path in args.comparison_tsvs:
        df = pd.read_csv(path, sep='\t')
        study = os.path.basename(path).split('_combined')[0]
        df['study'] = study
        dfs.append(df)

    all_df = pd.concat(dfs, ignore_index=True)
    all_df.to_csv('summary.classifier_comparison.tsv', sep='\t', index=False)

    # --- 1. Macro F1 bar chart ---
    macro_df = all_df[all_df['label'] == 'macro avg'].copy()
    mean_macro = (macro_df.groupby(['key', 'method'])['f1']
                  .agg(['mean', 'std']).reset_index())
    mean_macro.columns = ['key', 'method', 'mean_f1', 'std_f1']

    ordered_keys = [k for k in args.ref_keys if k in mean_macro['key'].unique()]
    pivot_mean = mean_macro.pivot(index='key', columns='method', values='mean_f1').loc[ordered_keys]
    pivot_std  = mean_macro.pivot(index='key', columns='method', values='std_f1').loc[ordered_keys]

    fig, ax = plt.subplots(figsize=(10, 7))
    x = np.arange(len(ordered_keys))
    width = 0.35
    methods = pivot_mean.columns.tolist()
    colors = {'RF': 'steelblue', 'kNN': 'coral'}
    for i, method in enumerate(methods):
        ax.bar(x + i * width, pivot_mean[method], width,
               yerr=pivot_std[method], capsize=4,
               label=method, color=colors.get(method, f'C{i}'), alpha=0.85)
    ax.set_xticks(x + width / 2)
    ax.set_xticklabels(ordered_keys, fontsize=13)
    ax.set_ylabel('Mean macro F1 (± std across studies)', fontsize=13)
    ax.set_title('RF vs kNN — mean macro F1 per hierarchy level', fontsize=14)
    ax.set_ylim(0, 1)
    ax.legend(title='Classifier', fontsize=12, title_fontsize=13)
    ax.tick_params(axis='y', labelsize=12)
    plt.tight_layout()
    plt.savefig('summary.classifier_comparison.macro_f1.png', dpi=150)
    plt.close()

    # --- 2. Per-class F1 heatmap per key ---
    per_class = all_df[~all_df['label'].isin(['macro avg', 'weighted avg', 'micro avg'])].copy()
    for key in args.ref_keys:
        key_df = per_class[per_class['key'] == key]
        if key_df.empty:
            continue

        mean_f1 = (key_df.groupby(['label', 'method'])['f1']
                   .mean().reset_index()
                   .pivot(index='label', columns='method', values='f1'))

        # order cell types by mean F1 across classifiers
        mean_f1 = mean_f1.loc[mean_f1.mean(axis=1).sort_values(ascending=False).index]

        fig_h = max(10, len(mean_f1) * 0.45)
        fig, ax = plt.subplots(figsize=(6, fig_h))
        sns.heatmap(mean_f1, ax=ax, vmin=0, vmax=1, cmap='YlOrRd',
                    xticklabels=True, yticklabels=True,
                    annot=True, fmt='.2f', annot_kws={'size': 11},
                    linewidths=0.4, linecolor='lightgrey')
        ax.set_title(f'Mean per-class F1 ({key})', fontsize=14)
        ax.set_xlabel('Classifier', fontsize=13)
        ax.set_ylabel('Cell type',  fontsize=13)
        plt.xticks(fontsize=13)
        plt.yticks(rotation=0, fontsize=11)
        ax.collections[0].colorbar.ax.tick_params(labelsize=11)
        plt.tight_layout()
        plt.savefig(f'summary.{key}.perclass_f1.png', dpi=150)
        plt.close()


if __name__ == '__main__':
    main()
