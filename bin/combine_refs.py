#!/usr/bin/env python3
"""
combine_refs.py

Concatenates all per-dataset ref h5ads from get_census_adata into a single
whole_cortex.h5ad for use as the combined reference in validation analyses.
"""

import argparse
import anndata as ad


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('--ref_paths', type=str, nargs='+', required=True)
    parser.add_argument('--output',    type=str, default='whole_cortex.h5ad')
    if __name__ == '__main__':
        known_args, _ = parser.parse_known_args()
        return known_args
    return parser.parse_args()


def main():
    args = parse_arguments()
    refs = [ad.read_h5ad(p) for p in args.ref_paths]
    combined = ad.concat(refs, merge='same')
    combined.obs_names_make_unique()
    combined.write_h5ad(args.output)


if __name__ == '__main__':
    main()
