#! /usr/bin/env Rscript


require(liana)
require(tibble)
require(purrr)
require(Seurat)
library(stringr)
library(dplyr)

sample_loc <- commandArgs(trailingOnly = TRUE)[1]
sample_id <- str_extract(str_extract(sample_loc, '[^/]+$'), '[^\\.]+')
print(paste('Working on', sample_id))


# Read in data
dat <- readRDS(sample_loc)

# Run liana
liana_dat <- liana_wrap(dat,
                         method = c('cellchat', 'squidpy', 'natmi', 'logfc', 'sca'),
                         resource = c('OmniPath'),
			 squidpy.params = list(cluster_key = 'ct'),
                         call_natmi.params = list(.natmi_path = '/data/sennis/AML/cc_interactions/liana/NATMI', num_cor = 8)
                         )


# Save results
saveRDS(liana_dat, file = str_replace(sample_loc, '.Rds', '_liana.Rds'))

# Run aggregate
dat_aggreg <- liana_dat %>% liana_aggregate()
saveRDS(dat_aggreg, file = str_replace(sample_loc, '.Rds', '_liana_aggreg.Rds'))
