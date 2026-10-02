#!/usr/bin/env python3
"""
get_combined_ref.py

Fetches reference data from CellXGene Census and writes the pre-built combined
reference (refs["whole cortex"]) directly to whole_cortex.h5ad.

get_census() already fetches all cells together as a single AnnData alongside
the per-dataset splits. This script just saves that combined object, avoiding
the need for a separate combine_refs step.
"""

import warnings
warnings.filterwarnings("ignore")

import os
import sys
import argparse
import utils
from utils import get_census


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('--organism',        type=str, required=True)
    parser.add_argument('--organ',           type=str, default='brain')
    parser.add_argument('--census_version',  type=str, required=True)
    parser.add_argument('--subsample_ref',   type=int, default=100)
    parser.add_argument('--relabel_path',    type=str, required=True)
    parser.add_argument('--ref_collections', type=str, nargs='+', required=True)
    parser.add_argument('--ref_keys',        type=str, nargs='+',
                        default=['subclass', 'class', 'family', 'global'])
    parser.add_argument('--seed',            type=int, default=42)
    if __name__ == '__main__':
        known_args, _ = parser.parse_known_args()
        return known_args
    return parser.parse_args()


def main():
    args = parse_arguments()

    if args.organism == 'mus_musculus':
        original_celltypes = utils.get_original_celltypes(
            columns_file=f"/space/grp/rschwartz/rschwartz/eval-references/meta/author_cell_annotations/{args.census_version}/original_celltype_columns.tsv",
            author_annotations_path=f"/space/grp/rschwartz/rschwartz/eval-references/meta/author_cell_annotations/{args.census_version}"
        )
    else:
        original_celltypes = None

    refs = get_census(
        organism=args.organism,
        organ=args.organ,
        subsample=args.subsample_ref,
        split_column='dataset_title',
        census_version=args.census_version,
        relabel_path=args.relabel_path,
        ref_collections=args.ref_collections,
        seed=args.seed,
        ref_keys=args.ref_keys,
        original_celltypes=original_celltypes
    )

    combined = refs.get("whole cortex")
    if combined is None:
        print("ERROR: 'whole cortex' key not found in refs dict.", file=sys.stderr)
        sys.exit(1)

    combined.write_h5ad('whole_cortex.h5ad')
    print(f"Written whole_cortex.h5ad ({combined.n_obs} cells, {combined.n_vars} genes)")


if __name__ == '__main__':
    main()
