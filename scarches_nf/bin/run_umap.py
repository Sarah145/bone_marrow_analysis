#!/usr/bin/env python3

# run scVI method

import argparse
import scanpy as sc
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# read in parameters
if __name__== "__main__":

    parser = argparse.ArgumentParser(description='Input/Output files and method parameters')

    parser.add_argument("--input_object",
                        dest='input_object',
                        type=str,
                        help ='Input h5ad file.')
    parser.add_argument("--output_prefix",
                        dest='output_prefix',
                        type=str,
                        help='Name to use as prefix for saving objects.')
    parser.add_argument("--batch_key",
                        dest='batch_key',
                        type=str,
                        default='Batch',
                        help ='Obs key defining batch.')
    parser.add_argument('--celltype_key',
                        dest='celltype_key',
                        type=str,
                        help='Obs key defining celltype.')
    parser.add_argument('--knn',
			dest='knn',
			type=int,
			help='Number of neighbours for computing knn graph.')
    args = parser.parse_args()

sc.settings.set_figure_params(dpi=200, frameon=False)
sc.set_figure_params(dpi=200)
sc.set_figure_params(figsize=(4, 4))
    
# read in data
adata = sc.read(args.input_object)

# compute knn and umap
sc.pp.neighbors(adata, n_neighbors=args.knn)
sc.tl.umap(adata)
sc.pl.umap(adata,
           color=[args.batch_key, args.celltype_key],
           frameon=False,
           wspace=0.6,
           ncols=1,
           save='_' + args.output_prefix + '.png',
           show=False)

df = pd.DataFrame(adata.obs)
umap = pd.DataFrame(list(adata.obsm['X_umap']), columns=['UMAP1', 'UMAP2'], index=adata.obs_names)
df['UMAP1'] = umap['UMAP1']
df['UMAP2'] = umap['UMAP2']

save_path = args.output_prefix + '_umap.csv'
df.to_csv(save_path, sep='\t')
