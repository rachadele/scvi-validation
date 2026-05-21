#!/usr/bin/env python3
"""
aggregate_query.py

Concatenates all processed query h5ads for a study and subsamples to n_cells.

Outputs:
  {study}_combined.h5ad
"""

import argparse
import os
import numpy as np
import anndata as ad


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('--query_paths', type=str, nargs='+', required=True)
    parser.add_argument('--study',       type=str, required=True)
    parser.add_argument('--n_cells',     type=int, default=1000)
    parser.add_argument('--seed',        type=int, default=42)
    if __name__ == '__main__':
        known_args, _ = parser.parse_known_args()
        return known_args
    return parser.parse_args()


def main():
    args = parse_arguments()
    np.random.seed(args.seed)

    adatas = [ad.read_h5ad(p) for p in args.query_paths]
    combined = ad.concat(adatas, label='sample', keys=[os.path.basename(p) for p in args.query_paths],
                         merge='same')

    if combined.n_obs > args.n_cells:
        idx = np.random.choice(combined.n_obs, size=args.n_cells, replace=False)
        combined = combined[idx, :]

    combined.obs.index = combined.obs.index.astype(str)
    combined.write_h5ad(f'{args.study}_combined.h5ad')


if __name__ == '__main__':
    main()
