library(tidyverse)
library(speckle)

# read in data and clean names
md <- read.csv('data/bone_marrow_integrated_umap.csv', sep = '\t', row.names = 1)
md <- md %>% rownames_to_column(var = 'cell') %>% 
  mutate(timepoint = case_when(
    timepoint == 'healthy' ~ 'Healthy',
    timepoint %in% c('DG', 'D0') ~ 'AML.Diagnosis',
    timepoint == 'R' ~ 'AML.Relapse',
    T ~ 'AML.Posttreatment'
  ),
  label = celltype_final)

# remove sorted samples
md <- md %>% filter(sample %ni% c('CD34_D2T1', 'CD34_D3T1', 'Normal_sorted_170531', 'Normal_sorted_170607', 'BM5-34p', 'BM5-34p38n'))

# set up list of comparisons
compars <- list(c('Healthy', 'AML.Diagnosis'), c('Healthy', 'AML.Relapse'), c('Healthy', 'AML.Posttreatment'),
                c('AML.Diagnosis', 'AML.Posttreatment'), c('AML.Diagnosis', 'AML.Relapse'),
                c('AML.Posttreatment', 'AML.Relapse'))

# run propeller for each comparison
res <- data.frame()
for(i in compars){
  sub <- md %>% filter(timepoint %in% i)
  prop <- propeller(clusters = sub$label, sample = sub$sample, 
                    group = sub$timepoint)
  prop$compar <- paste0(i, collapse = ', ')
  prop <- prop %>% pivot_longer(starts_with('PropMean'), names_to = 'tp', values_to = 'prop') %>% mutate(tp = str_sub(tp, start = 10))
  res <- rbind(res, prop)
}

p_df <- res %>% separate(compar, into = c('TP1', 'TP2'), sep = ', ')

TP1_greater <- c()
for(i in 1:nrow(p_df)){
  sub <- p_df %>% filter(TP1 == p_df$TP1[i], TP2 == p_df$TP2[i], BaselineProp.clusters == p_df$BaselineProp.clusters[i])
  TP1_greater <- c(TP1_greater, (sub %>% filter(tp == TP1) %>% pull(prop)) > (sub %>% filter(tp == TP2) %>% pull(prop)))
}
p_df$TP1_greater <- TP1_greater


# take quick look at results
ggplot(p_df %>% select(tp, prop, ct = BaselineProp.clusters) %>% mutate(prop = round(prop, 5)) %>% distinct(),
       aes(x = tp, y = prop, fill = ct)) +
  geom_col()

ggplot(p_df %>% filter(FDR <= 0.05), aes(x = TP2, y = TP1, size = -log10(FDR), col = TP1_greater)) + 
  geom_point() + 
  facet_wrap(~BaselineProp.clusters)

# save results
saveRDS(p_df, 'data/propeller_results.Rds')
