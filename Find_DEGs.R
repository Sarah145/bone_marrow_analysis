# Load libraries
library(tidyverse)
library(Seurat)
library(Libra)
library(clusterProfiler)
library(org.Hs.eg.db)

# Read in data (subset containing only healthy and AML samples)
dat <- readRDS('../data/healthy_d0_seu_object.Rds')

# Identify DEGs using libra
degs <- run_de(dat, cell_type_col = 'celltype_final', replicate_col = 'sample', label_col = 'tp', min_features = 10)
degs$sig <- ifelse(abs(degs$avg_logFC) >= 1 & degs$p_val_adj <= 0.05, T, F)
saveRDS(degs, '../data/d0_healthy_degs_libra.Rds')

# Get background gene list for each cell type
bg_gene_list <- list()
for(i in unique(degs$cell_type)){
  sub <- dat[,dat$celltype_final == i]
  bg_genes <- names(rowSums(sub@assays$RNA@counts)[rowSums(sub@assays$RNA@counts) > 0])
  bg_gene_list[[i]] <- bg_genes
}

# Identify enriched GO terms for each cell type
d0_go_df <- data.frame()
healthy_go_df <- data.frame()
for(i in unique(degs$cell_type)){
  up_de_genes <- degs %>% filter(cell_type == i, sig == T, avg_logFC > 0) %>% pull(gene)
  down_de_genes <- degs %>% filter(cell_type == i, sig == T, avg_logFC < 0) %>% pull(gene)
  up_go <- data.frame(enrichGO(gene = up_de_genes, universe = bg_gene_list[[i]], OrgDb = org.Hs.eg.db, keyType = 'SYMBOL', ont = "BP")) %>% mutate(cell_type = i, timepoint = 'AML (diagnosis)')
  down_go <- data.frame(enrichGO(gene = down_de_genes, universe = bg_gene_list[[i]], OrgDb = org.Hs.eg.db, keyType = 'SYMBOL', ont = "BP")) %>% mutate(cell_type = i, timepoint = 'Healthy')
  d0_go_df <- rbind(d0_go_df, up_go)
  healthy_go_df <- rbind(healthy_go_df, down_go)
}

all_gos <- rbind(d0_go_df, healthy_go_df)
saveRDS(all_gos, '../data/healthy_d0_enriched_gos_1FC.Rds')
