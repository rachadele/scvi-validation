#!/usr/bin/env python3
"""
compare_classifiers.py

Trains RF (subclass only, via utils.rfc_pred) and weighted kNN on ref obsm["scvi"],
predicts subclass on query, then maps predictions up the hierarchy using aggregate_labels.
Computes per-class and macro F1 at every hierarchy level for both classifiers.

PYTHONPATH must include bin/ so utils can be imported (handled by symlink).

Outputs:
  {query_name}_{ref_name}.classifier_comparison.tsv
  {query_name}_{ref_name}.classifier_comparison.png
"""

import argparse
import os
import sys
import numpy as np
import pandas as pd
import anndata as ad
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from sklearn.ensemble import RandomForestClassifier
from sklearn.neighbors import KNeighborsClassifier
from sklearn.metrics import classification_report
import random


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('--ref_path',     type=str, required=True)
    parser.add_argument('--query_path',   type=str, required=True)
    parser.add_argument('--mapping_file', type=str, required=True)
    parser.add_argument('--ref_keys',     type=str, nargs='+',
                        default=['subclass', 'class', 'family', 'global'])
    parser.add_argument('--n_neighbors',  type=int, default=15)
    parser.add_argument('--seed',         type=int, default=42)
    if __name__ == '__main__':
        known_args, _ = parser.parse_known_args()
        return known_args
    return parser.parse_args()


def collect_f1(y_true, y_pred, method, key):
    # Only evaluate labels present in the query (support > 0)
    present_labels = sorted(set(y_true))
    report = classification_report(y_true, y_pred, labels=present_labels,
                                   output_dict=True, zero_division=np.nan)
    records = []
    for label, metrics in report.items():
        if label in ('accuracy', 'macro avg', 'weighted avg'):
            continue
        records.append({
            'method':  method,
            'key':     key,
            'label':   label,
            'f1':      metrics['f1-score'],
            'support': metrics['support'],
        })
    records.append({
        'method':  method,
        'key':     key,
        'label':   'macro avg',
        'f1':      report['macro avg']['f1-score'],
        'support': report['macro avg']['support'],
    })
    return records


def main():
    args = parse_arguments()
    random.seed(args.seed)
    np.random.seed(args.seed)

    from utils import aggregate_labels, map_valid_labels
    from sklearn.ensemble import RandomForestClassifier

    ref_name   = os.path.basename(args.ref_path).replace('.h5ad', '')
    query_name = os.path.basename(args.query_path).replace('.h5ad', '')
    prefix     = f'{query_name}_{ref_name}'

    ref   = ad.read_h5ad(args.ref_path,   backed='r')
    query = ad.read_h5ad(args.query_path, backed='r')
    mapping_df = pd.read_csv(args.mapping_file, sep='\t')

    granular_key = args.ref_keys[0]
    X_ref   = ref.obsm['scvi']
    X_query = query.obsm['scvi']
    y_ref   = ref.obs[granular_key].astype(str).values

    # Base DataFrame of true labels from query
    query_obs = query.obs.reset_index(drop=True)

    def build_pred_df(y_subclass_pred):
        """Aggregate predictions up hierarchy then fix coarse query labels via map_valid_labels."""
        pred_df = pd.DataFrame({'predicted_' + granular_key: y_subclass_pred},
                               index=range(len(y_subclass_pred)))
        pred_df = aggregate_labels(pred_df, mapping_df, args.ref_keys, predicted=True)
        # Combine with true labels so map_valid_labels can promote coarse query labels
        combined = query_obs[[k for k in args.ref_keys if k in query_obs.columns]].copy()
        combined = combined.reset_index(drop=True)
        for col in pred_df.columns:
            combined[col] = pred_df[col].values
        combined = map_valid_labels(combined, args.ref_keys, mapping_df)
        return combined

    # --- RF: train at subclass, aggregate up hierarchy ---
    rfc = RandomForestClassifier(class_weight='balanced', random_state=args.seed)
    rfc.fit(X_ref, y_ref)
    prob_mat      = rfc.predict_proba(X_query)
    class_labels  = np.array(rfc.classes_)
    y_rf_subclass = class_labels[np.argmax(prob_mat, axis=1)]
    rf_combined   = build_pred_df(y_rf_subclass)

    # --- kNN: train at subclass, aggregate up hierarchy ---
    knn = KNeighborsClassifier(n_neighbors=args.n_neighbors, weights='distance',
                               metric='cosine', n_jobs=-1)
    knn.fit(X_ref, y_ref)
    y_knn_subclass = knn.predict(X_query)
    knn_combined   = build_pred_df(y_knn_subclass)

    # --- Compute F1 at every hierarchy level ---
    records = []
    for key in args.ref_keys:
        if key not in rf_combined.columns:
            continue
        y_true   = rf_combined[key].astype(str).values
        pred_col = 'predicted_' + key

        if pred_col in rf_combined.columns:
            records.extend(collect_f1(y_true, rf_combined[pred_col].astype(str).values, 'RF', key))
        if pred_col in knn_combined.columns:
            records.extend(collect_f1(y_true, knn_combined[pred_col].astype(str).values, 'kNN', key))

    df = pd.DataFrame(records)
    df.to_csv(f'{prefix}.classifier_comparison.tsv', sep='\t', index=False)

    # Bar chart: macro F1 per hierarchy level for RF vs kNN
    macro_df = df[df['label'] == 'macro avg'].copy()
    if not macro_df.empty:
        pivot = macro_df.pivot(index='key', columns='method', values='f1')
        # preserve hierarchy order
        ordered = [k for k in args.ref_keys if k in pivot.index]
        pivot = pivot.loc[ordered]
        fig, ax = plt.subplots(figsize=(10, 7))
        pivot.plot(kind='bar', ax=ax, rot=0, color=['steelblue', 'coral'])
        ax.set_ylabel('Macro F1')
        ax.set_title(f'RF vs kNN — {query_name} × {ref_name}', fontsize=8)
        ax.set_ylim(0, 1)
        ax.legend(title='Classifier')
        plt.tight_layout()
        plt.savefig(f'{prefix}.classifier_comparison.png', dpi=150)
        plt.close()


if __name__ == '__main__':
    main()
