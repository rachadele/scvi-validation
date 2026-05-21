#!/usr/bin/env python3
"""
joint_umap.py

Concatenates ref and query on obsm["scvi"], computes a joint UMAP, and saves plots:
  - colored by cell type (ref_keys[0])
  - colored by origin (ref vs query)

Outputs:
  {query_name}_{ref_name}.celltype.umap.png
  {query_name}_{ref_name}.origin.umap.png
"""

import argparse
import os
import sys
import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('--ref_path',   type=str, required=True)
    parser.add_argument('--query_path', type=str, required=True)
    parser.add_argument('--ref_keys',   type=str, nargs='+',
                        default=['subclass', 'class', 'family', 'global'])
    parser.add_argument('--n_neighbors', type=int, default=15)
    parser.add_argument('--seed', type=int, default=42)
    if __name__ == '__main__':
        known_args, _ = parser.parse_known_args()
        return known_args
    return parser.parse_args()


def main():
    args = parse_arguments()
    ref_name   = os.path.basename(args.ref_path).replace('.h5ad', '')
    query_name = os.path.basename(args.query_path).replace('_processed.h5ad', '')
    prefix     = f'{query_name}_{ref_name}'

    ref   = ad.read_h5ad(args.ref_path)
    query = ad.read_h5ad(args.query_path)

    # Keep only obs columns present in both
    shared_cols = list(set(ref.obs.columns) & set(query.obs.columns))
    ref_sub   = ad.AnnData(X=ref.X,   obs=ref.obs[shared_cols],   obsm={'scvi': ref.obsm['scvi']})
    query_sub = ad.AnnData(X=query.X, obs=query.obs[shared_cols], obsm={'scvi': query.obsm['scvi']})

    ref_sub.obs['_origin']   = 'ref'
    query_sub.obs['_origin'] = 'query'

    combined = ad.concat([ref_sub, query_sub], label='_source', keys=['ref', 'query'])
    combined.obsm['scvi'] = np.vstack([ref.obsm['scvi'], query.obsm['scvi']])

    sc.pp.neighbors(combined, use_rep='scvi', n_neighbors=args.n_neighbors,
                    random_state=args.seed)
    sc.tl.umap(combined, random_state=args.seed)

    # Plot by cell type (ref_keys[0])
    ct_key = args.ref_keys[0]
    if ct_key in combined.obs.columns:
        plt.figure(figsize=(12, 10))
        fig = sc.pl.umap(combined, color=ct_key, title=f'{prefix} — {ct_key}',
                         show=False, return_fig=True, legend_loc='on data',
                         legend_fontsize=6)
        fig.set_size_inches(12, 10)
        fig.savefig(f'{prefix}.celltype.umap.png', dpi=150, bbox_inches='tight')
        plt.close(fig)

    # Plot by origin
    plt.figure(figsize=(12, 10))
    fig = sc.pl.umap(combined, color='_origin', title=f'{prefix} — origin',
                     show=False, return_fig=True)
    fig.set_size_inches(12, 10)
    fig.savefig(f'{prefix}.origin.umap.png', dpi=150, bbox_inches='tight')
    plt.close(fig)


if __name__ == '__main__':
    main()
