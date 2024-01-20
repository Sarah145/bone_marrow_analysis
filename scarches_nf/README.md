## Nextflow pipeline for batch correction

This pipeline runs scArches (scVI) to integrate multiple single-cell datasets, trying out different values for certain hyperparameters (number of latent dimensions, number of highly variable genes and number of nearest neighbours for the UMAP graph), and then evaluating each output using the scPOP R package to find the optimal hyperparameters.
