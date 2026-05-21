#!/usr/bin/env python3
"""
plot_embedding_coherence.py

Reads centroid cosine similarity TSV files from both scVI and Seurat coherence
steps and produces heatmap PNGs organised by method and study:
  {method}/{study}/{prefix}.{key}.{suffix}.png

Handles both file types:
  {prefix}.{key}.coherence.tsv        -> method=scvi
  {prefix}.{key}.seurat_coherence.tsv -> method=seurat
"""

import argparse
import os
import re
import sys
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from coherence_utils import plot_heatmap


def parse_arguments():
    parser = argparse.ArgumentParser()
    parser.add_argument('--coherence_tsvs', type=str, nargs='+', required=True)
    if __name__ == '__main__':
        known_args, _ = parser.parse_known_args()
        return known_args
    return parser.parse_args()


def parse_filename(path):
    """Return (prefix, study, key, method, out_png) from a coherence TSV filename."""
    name = os.path.basename(path)
    m = re.match(r'^(.+)\.(\w+)\.(seurat_coherence|coherence)\.tsv$', name)
    if not m:
        print(f'WARNING: could not parse filename {name}, skipping', file=sys.stderr)
        return None
    prefix, key, suffix = m.group(1), m.group(2), m.group(3)
    method  = 'seurat' if suffix == 'seurat_coherence' else 'scvi'
    study   = prefix.split('_combined')[0]
    out_png = re.sub(r'\.tsv$', '.png', name)
    return prefix, study, key, method, out_png


def main():
    args = parse_arguments()

    for tsv_path in args.coherence_tsvs:
        parsed = parse_filename(tsv_path)
        if parsed is None:
            continue
        prefix, study, key, method, out_png = parsed

        parts      = prefix.split('_', 1)
        query_name = parts[0] if len(parts) == 2 else prefix
        ref_name   = parts[1] if len(parts) == 2 else ''

        out_dir = os.path.join(method, study)
        os.makedirs(out_dir, exist_ok=True)
        out_path = os.path.join(out_dir, out_png)

        df = pd.read_csv(tsv_path, sep='\t', index_col=0)
        df.index.name   = 'ref_cell_type'
        df.columns.name = 'query_cell_type'

        title = f'Centroid cosine similarity [{method.upper()}] ({key})\n{query_name} vs {ref_name}'
        plot_heatmap(df, title, out_path)
        print(f'Wrote {out_path}')


if __name__ == '__main__':
    main()
