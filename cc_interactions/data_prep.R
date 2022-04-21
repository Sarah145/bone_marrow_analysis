# Load libraries
library(SeuratDisk)
library(Seurat)
library(tidyverse)

# need to convert h5ad files to seurat objects to run analysis
Convert('../data/healthy_sub.h5ad', dest = 'h5seurat', overwrite = F)
Convert('../data/d0_sub.h5ad', dest = 'h5seurat', overwrite = F)

hea <- LoadH5Seurat('../data/healthy_sub.h5seurat', assay = 'RNA')
d0 <- LoadH5Seurat('../data/d0_sub.h5seurat', assay = 'RNA')

# need to run FindVariableFeatures or else squidpy won't run
hea <- FindVariableFeatures(hea, nfeatures = nrow(hea))
d0 <- FindVariableFeatures(d0, nfeatures = nrow(d0))

# cell type names can't have dots in them
hea$ct <- str_remove(hea$celltype_final, '\\.')
d0$ct <- str_remove(d0$celltype_final, '\\.')

# convert names to factor and set as Idents
hea$ct <- factor(hea$ct)
d0$ct <- factor(d0$ct)
Idents(hea) <- hea$ct
Idents(d0) <- d0$ct

# save files to run interaction prediction
saveRDS(hea, '../data/healthy_sub.Rds')
saveRDS(d0, '../data/d0_sub.Rds')
