#!/usr/bin/env python3
"""
process_query_val.py — like process_query.py but with per-dataset subsampling support.

Extra args vs the original:
  --subsample_per_dataset N   : sample N cells per unique value of --dataset_key (default: sample_id)
  --dataset_key COL           : obs column to group by for per-dataset subsampling
  --eval_pipeline_dir PATH    : path to nextflow_eval_pipeline (to import utils.py)
"""

import argparse
import os
import sys
import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
import scvi


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('--model_path', type=str, required=True)
    parser.add_argument('--relabel_path', type=str, required=True)
    parser.add_argument('--query_path', type=str, required=True)
    parser.add_argument('--mapping_file', type=str, default=None)
    parser.add_argument('--batch_key', type=str, default='sample_id')
    parser.add_argument('--join_key', type=str, default=None)
    parser.add_argument('--ref_keys', type=str, nargs='+', default=['subclass', 'class', 'family', 'global'])
    parser.add_argument('--remove_unknown', action='store_true')
    parser.add_argument('--nmads', type=int, default=5)
    parser.add_argument('--gene_mapping', type=str, required=True)
    parser.add_argument('--seed', type=int, default=42)
    # subsampling — mutually exclusive approaches
    parser.add_argument('--subsample_query', type=int, default=None,
                        help='Global subsample: keep this many cells per query file')
    parser.add_argument('--subsample_per_dataset', type=int, default=None,
                        help='Per-dataset subsample: keep this many cells per unique value of --dataset_key')
    parser.add_argument('--dataset_key', type=str, default='study',
                        help='obs column to group by for per-dataset subsampling')
    if __name__ == '__main__':
        known_args, _ = parser.parse_known_args()
        return known_args
    return parser.parse_args()


def subsample_global(query, n, seed):
    if query.n_obs > n:
        idx = np.random.choice(query.n_obs, size=n, replace=False)
        query = query[idx, :]
        query.obs.index = query.obs.index.astype(str)
    return query


def subsample_per_dataset(query, n, dataset_key, seed):
    """Sample up to n cells per unique value of dataset_key in obs."""
    if dataset_key not in query.obs.columns:
        print(f"Warning: '{dataset_key}' not in obs columns, falling back to global subsample of {n}",
              file=sys.stderr)
        return subsample_global(query, n, seed)

    rng = np.random.default_rng(seed)
    keep = []
    for _, grp in query.obs.groupby(dataset_key):
        n_take = min(len(grp), n)
        keep.extend(rng.choice(grp.index, size=n_take, replace=False).tolist())
    query = query[keep, :]
    query.obs.index = query.obs.index.astype(str)
    return query


def main():
    args = parse_arguments()

    import utils
    from utils import map_genes, relabel, aggregate_labels, get_qc_metrics, process_query
    import random

    SEED = args.seed
    random.seed(SEED)
    np.random.seed(SEED)
    scvi.settings.seed = SEED

    query = ad.read_h5ad(args.query_path)

    # Subsampling — per-dataset takes priority over global
    if args.subsample_per_dataset is not None:
        query = subsample_per_dataset(query, args.subsample_per_dataset, args.dataset_key, SEED)
    elif args.subsample_query is not None:
        query = subsample_global(query, args.subsample_query, SEED)

    sc.pp.scrublet(query, batch_key=None)

    gene_mapping = pd.read_csv(args.gene_mapping, sep='\t', header=0)
    query = map_genes(query, gene_mapping)

    mapping_df = pd.read_csv(args.mapping_file, sep='\t')
    query = relabel(query=query, relabel_path=args.relabel_path, join_key=args.join_key, sep='\t')

    nan_mask = query.obs['subclass'].isna()
    if nan_mask.any():
        raise ValueError(
            f"relabeling produced NaN in 'subclass' for {nan_mask.sum()} cells. "
            f"Check that all cell type labels in the query are covered by {args.relabel_path}."
        )

    query.obs = aggregate_labels(query=query.obs, mapping_df=mapping_df,
                                 ref_keys=args.ref_keys, predicted=False)

    if args.remove_unknown:
        query = query[query.obs[args.ref_keys[0]] != 'unknown']

    query.obs.index = range(query.n_obs)
    query = get_qc_metrics(query, args.nmads)

    raw_name = os.path.basename(args.query_path).replace('.h5ad', '_raw')
    for col in query.obs.select_dtypes(include=['object', 'category']).columns:
        query.obs[col] = query.obs[col].astype(str)
    query.write_h5ad(f'{raw_name}.h5ad')

    query = process_query(query, args.model_path, args.batch_key, seed=SEED)

    processed_name = os.path.basename(args.query_path).replace('.h5ad', '_processed')
    query.write_h5ad(f'{processed_name}.h5ad')


if __name__ == '__main__':
    main()
