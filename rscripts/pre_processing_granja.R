# Processing granja 2019 data

# ---------

library(tidyr)
library(dplyr)
library(ggplot2)
library(tibble)
library(Matrix)
library(stringr)

# Start with scADT data
# raw data available as rds file from GEO under accession GSE139369
# only interested in healthy data so just downloaded GSM4138872, GSM4138873, GSM4138874, GSM4138875	
# download link for processed data is also available on their GitHub but that rds file only contains 6 proteins?

# read in data from GEO
f1 <- t(readRDS('../data/granja/GSM4138872_scADT_BMMC_D1T1.rds'))
f2 <- t(readRDS('../data/granja/GSM4138873_scADT_BMMC_D1T2.rds'))
f3 <- t(readRDS('../data/granja/GSM4138874_scADT_CD34_D2T1.rds'))
f4 <- t(readRDS('../data/granja/GSM4138875_scADT_CD34_D3T1.rds'))
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
write.table(gran_clr, '../data/granja/gran_adt_clr.csv', sep = '\t', quote = F)

#-------

# Now move on the scRNA data
# The rds that was available from their github was not right - possibly mislabelled cells? - so went with raw counts available from GEO
# read in data from GEO

f1 <- readRDS('../data/granja/GSM4138872_scRNA_BMMC_D1T1.rds')
f2 <- readRDS('../data/granja/GSM4138873_scRNA_BMMC_D1T2.rds')
f3 <- readRDS('../data/granja/GSM4138874_scRNA_CD34_D2T1.rds')
f4 <- readRDS('../data/granja/GSM4138875_scRNA_CD34_D3T1.rds')
gran_rna <- cbind(f1, f2, f3, f4)
dim(gran_rna)

# need to write out 3 separate files so I can read it in to python

# write out counts matrix
writeMM(t(gran_rna), '../data/granja/gran_rna.mtx')

# write out gene names
gns <- data.frame(gene_symbol = rownames(gran_rna))
write.table(gns, '../data/granja/gran_rna_genes.csv', sep = '\t', quote = F, row.names = F)

# write out barcode names and metadata
# rds file on their GitHub has metadata for all cells - download link for the rds file here: https://jeffgranja.s3.amazonaws.com/MPAL-10x/Supplementary_Data/Healthy-Data/scRNA-Healthy-Hematopoiesis-191120.rds
processed_gran <- readRDS('../data/granja/scRNA-Healthy-Hematopoiesis-191120.rds')
meta <- data.frame(processed_gran@colData)
meta$bc <- paste0(meta$Group, ':', meta$Barcode)
rownames(meta) <- meta$bc
meta <- meta[colnames(gran_rna), ]
write.table(meta, '../data/granja/gran_rna_meta.csv', sep = '\t', quote = F)

# I'll do the rest of the pre-processing in python ('../notebooks/pre_processing_granja.ipynb')