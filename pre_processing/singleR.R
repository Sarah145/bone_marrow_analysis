#! /usr/bin/env Rscript


library(SingleCellExperiment)
library(SingleR)
library(stringr)

ref_loc <- commandArgs(trailingOnly = TRUE)[1]
sample_loc <- commandArgs(trailingOnly = TRUE)[2]
sample_id <- str_extract(str_extract(sample_loc, '[^/]+$'), '[^\\.]+')
print(paste('Working on', sample_id))

# set up function to convert h5ad to sce
#function to build sce
build_sce <- function(counts_mat, meta_data, assay_name, emb_list){
  sce <- SingleCellExperiment(assays = list(RNA = counts_mat),
                              colData = meta_data,
                              rowData = list('feature_names' = rownames(counts_mat)),
                              reducedDims = emb_list)
  sce
}



h5ad2sce <- function(h5ad_path){
  anndata <- reticulate::import('anndata', convert = F)
  pandas <- reticulate::import('pandas', convert = F)
  numpy <- reticulate::import('numpy', convert = F) 
  #Read h5ad file
  h5ad_file <- anndata$read_h5ad(h5ad_path)
  #Extract Counts   
  #set up condition ERROR this is because logcounts and scanorama counts are numpy objects with different properties
  h5ad_counts <- tryCatch({
    pandas$DataFrame(h5ad_file$X, index = h5ad_file$obs_names, columns = h5ad_file$var_names)
  }, error=function(e){pandas$DataFrame(h5ad_file$X$todense(), index = h5ad_file$obs_names, columns = h5ad_file$var_names)})
  
  counts_mat <- t(as.matrix(reticulate::py_to_r(h5ad_counts)))
  #Extract Metadata
  meta_data <- reticulate::py_to_r(pandas$DataFrame(h5ad_file$obs, dtype = "object"))
  #Extract Embeddings
  emb_names <- reticulate::py_to_r(h5ad_file$obsm_keys())
  emb_list <- lapply(as.list(emb_names), function(x) reticulate::py_to_r(h5ad_file$obsm[x]))
  names(emb_list) <- emb_names
  
  #Build SCE
  sce <- build_sce(counts_mat = counts_mat, meta_data = meta_data, emb_list = emb_list)
  sce
}


# read in data as sce to use as reference
ref <- h5ad2sce(ref_loc)

# read in sample data as sce
sam <- h5ad2sce(sample_loc)

# run SingleR
pred <- SingleR(test=sam, ref=ref, labels=ref[['celltype']], de.method="wilcox", assay.type.test = 'RNA', assay.type.ref = 'RNA')

df <- data.frame(label = pred$labels)
rownames(df) <- rownames(pred)

# save results
saveRDS(pred, file = str_replace(sample_loc, '.h5ad', '_singleR_results.rds'))
write.table(df, file = str_replace(sample_loc, '.h5ad', '_cell_labels.csv'), sep = '\t', row.names = T, col.names = T, quote = F)

print('Done!')
