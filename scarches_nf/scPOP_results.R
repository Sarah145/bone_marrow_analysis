library(tidyverse)
library(scPOP)

tbl <-
  list.files(pattern = "data/eval/*.csv") %>% 
  map_df(~read.csv(., sep = '\t')) 
tbl <- tbl %>% select(-lisi_celltype.1, -silWidth_celltype.1)

# there's a bug in the calc_sumZscore function so I copied the source code
# from github (https://github.com/vinay-swamy/scPOP/blob/main/R/metrics.R) 
# and fixed bug at line 274

calc_sumZscore <- function(df, batch_key){
  id <- df[,1]
  df <- df[,-1]
  batch_cols <- colnames(df)[grepl(batch_key, colnames(df))]
  label_cols <- colnames(df)[!grepl(batch_key, colnames(df))]
  negative_cases <- c(  batch_cols[grepl('silWidth', batch_cols)],
                        label_cols[grepl('lisi', label_cols)]
  )
  positive_cases <- c(  batch_cols[!grepl('silWidth', batch_cols)],
                        label_cols[!grepl('lisi', label_cols)]
  )
  nc_scaled <- apply(df[,negative_cases], 2, scale) * -1
  pc_scaled <- apply(df[,positive_cases], 2, scale)
  pc_scaled[is.nan(pc_scaled)] <- 0
  nc_scaled[is.nan(nc_scaled)] <- 0
  sumZ <- rowSums(cbind(pc_scaled, nc_scaled))
  return(sumZ)
}

tbl$sumZscore <-  calc_sumZscore(tbl, 'sample')

get_scaled_mtx <- function(df, batch_key){
  id <- df[,1]
  df <- df[,-1]
  batch_cols <- colnames(df)[grepl(batch_key, colnames(df))]
  label_cols <- colnames(df)[!grepl(batch_key, colnames(df))]
  negative_cases <- c(  batch_cols[grepl('silWidth', batch_cols)],
                        label_cols[grepl('lisi', label_cols)]
  )
  positive_cases <- c(  batch_cols[!grepl('silWidth', batch_cols)],
                        label_cols[!grepl('lisi', label_cols)]
  )
  nc_scaled <- apply(df[,negative_cases], 2, scale) * -1
  pc_scaled <- apply(df[,positive_cases], 2, scale)
  pc_scaled[is.nan(pc_scaled)] <- 0
  nc_scaled[is.nan(nc_scaled)] <- 0
  scaled_mtx <- cbind(pc_scaled, nc_scaled)
  return(scaled_mtx)
}

sc_mtx <- get_scaled_mtx(tbl %>% select(-sumZscore), 'sample') 

sc_tbl <- data.frame(cbind(run = tbl$run, sumZscore = as.numeric(tbl$sumZscore), sc_mtx)) %>% mutate(run = fct_reorder(run, as.numeric(sumZscore))) %>% pivot_longer(2:8, names_to = 'metric')

ggplot(sc_tbl, aes(x = metric, y = run, fill = as.numeric(value))) +
  geom_tile() +
  scale_x_discrete(limits = c('ari_label', 'nmi_label', 'lisi_sample', 'lisi_celltype', 'silWidth_sample', 'silWidth_celltype', 'sumZscore')) +
  scale_fill_viridis_c(option = 'H', direction = 1) +
  labs(x = 'Metric', y = 'Run', fill = 'Value', title = 'scPOP Results') +
  theme_classic(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text = element_text(colour = 'black'),
        axis.line = element_blank())
ggsave('scPOP_results.png', height = 12, width = 7)

write.table(tbl, 'scPOP_results.csv', sep='\t', row.names = F, quote = F)
