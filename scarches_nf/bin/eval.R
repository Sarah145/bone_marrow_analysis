#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))

option_list = list(
  make_option(
    c("-i", "--input_object"),
    action = "store",
    default = NA,
    type = 'character',
    help = 'Path to input file.'
  ),
  make_option(
    c("-o", "--output_prefix"),
    action = "store",
    default = NA,
    type = 'character',
    help = 'Prefix for naming output file.'
  ),
  make_option(
    c("-b", "--batch_key"),
    action = "store",
    default = "batch",
    type = 'character',
    help = 'Column defining batch.'
  ),
  make_option(
    c("-c", "--celltype_key"),
    action = "store",
    default = "celltype",
    type = 'character',
    help = 'Column defining celltype.'
  )
)

opt <- parse_args(OptionParser(option_list=option_list))

library(scPOP)


df <- read.csv(opt$input_object, sep = '\t', row.names = 1)
metrics <- run_all_metrics(reduction = df[,c('UMAP1', 'UMAP2')],
                           metadata = df,
                           batch_key = opt$batch_key,
                           label1_key = opt$celltype_key,
                           label2_key = opt$celltype_key,
                           sil_width_prop = 0.25,
                           run_name = opt$output_prefix)

write.table(metrics, file = paste0(opt$output_prefix, '_eval.csv'), sep='\t', row.names=F, quote=F)

