# Processing granja 2019 data

# ---------

library(tidyr)
library(dplyr)
library(ggplot2)
library(tibble)
library(here)
library(Matrix)
library(stringr)

# Start with scADT data
# raw data available as rds file from GEO under accession GSE139369
# only interested in healthy data so just downloaded GSM4138872, GSM4138873, GSM4138874, GSM4138875	
# download link processed data is also available on their GitHub but that rds file only contains 6 proteins?

# read in data from GEO
f1 <- t(readRDS(here('../data/granja/GSM4138872_scADT_BMMC_D1T1.rds')))
f2 <- t(readRDS(here('../data/granja/GSM4138873_scADT_BMMC_D1T2.rds')))
f3 <- t(readRDS(here('../data/granja/GSM4138874_scADT_CD34_D2T1.rds')))
f4 <- t(readRDS(here('../data/granja/GSM4138875_scADT_CD34_D3T1.rds')))
gran_adt <- data.frame(rbind(f1, f2, f3, f4))
dim(gran_adt)

# these are raw counts, need to convert to CLR-transformed values
# for clr function I adapted code from Seurat GitHub https://github.com/satijalab/seurat/blob/b56d194939379460db23380426d3896b54d91ab6/R/preprocessing.R

CLR <- function(x){
  return(log1p(x = x / (exp(x = sum(log1p(x = x[x > 0]), na.rm = TRUE) / length(x = x)))))
}

gran_clr <- t(apply(gran_adt, 1, CLR))

# update colnames of matrix to avoid confusion between protein and genes
colnames(gran_clr) <- paste0('ADT_', colnames(gran_clr))

# take quick look at the distribution of expression values for each protein
pd <- data.frame(gran_clr) %>% rownames_to_column('cell') %>% pivot_longer(1:ncol(gran_clr) +1, names_to = 'protein')
pd[is.na(pd)] <- 0
ggplot(pd, aes(x = value)) +
  geom_histogram(col = 'gray20', fill = 'deeppink', alpha = 0.6, bins=40) +
  facet_wrap(~protein, scales = 'free', nrow = 7) +
  theme_minimal(base_size = 14)

# write out matrix of clr-transformed adt counts
write.table(gran_clr, here('../data/granja/gran_adt_clr.csv'), sep = '\t', quote = F)

#-------

# Now move on the scRNA data
# I chose to use the rds that was available from their github because it has less cells which presumably means it's been QC-ed
# download link for the rds file here: https://jeffgranja.s3.amazonaws.com/MPAL-10x/Supplementary_Data/Healthy-Data/scRNA-Healthy-Hematopoiesis-191120.rds

gran_rna <- readRDS(here('../data/granja/scRNA-Healthy-Hematopoiesis-191120.rds'))

# I only want bone marrow so will exclude PBMC samples
gran_rna <- gran_rna[,!str_detect(colnames(gran_rna), 'PBMC')]

# need to write out 3 separate files so I can read it in to python

# write out counts matrix
mat <- t(gran_rna@assays$data$counts)
writeMM(mat, here('../data/granja/gran_rna.mtx'))

# write out gene names
gns <- data.frame(gene_symbol = colnames(mat))
write.table(gns, here('../data/granja/gran_rna_genes.csv'), sep = '\t', quote = F, row.names = F)

# write out barcode names and metadata
meta <- data.frame(gran_rna@colData)
write.table(meta, here('../data/granja/gran_rna_meta.csv'), sep = '\t', quote = F)

# I'll do the rest of the pre-processing in python ('../notebooks/pre_processing_granja.ipynb')