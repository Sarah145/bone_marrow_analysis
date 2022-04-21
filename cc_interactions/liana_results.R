library(tidyverse)
library(AnnotationDbi)
#library(org.Hs.eg.db)
library(OmnipathR)

# read in liana results and filter to only interactions with rank <= 0.05
# filter to only interactions involving HSC.MPPs
healthy <- readRDS('healthy_sub_liana_aggreg.Rds')
healthy_sig <- healthy %>% filter(aggregate_rank <= 0.05, (source %in% c('HSCMPP') | target %in% c('HSCMPP'))) %>% 
  mutate(tp = 'Healthy')

d0 <- readRDS('d0_sub_liana_aggreg.Rds')
d0_sig <- d0 %>% filter(aggregate_rank <= 0.05, (source %in% c('HSCMPP') | target %in% c('HSCMPP'))) %>% 
  mutate(tp = 'Diagnosis')

# combine results for both timepoints and fix cell type names
both_sig <- rbind(healthy_sig, d0_sig) %>%
  mutate(source = case_when(
    source == 'HSCMPP' ~ 'HSC.MPP',
    str_detect(source, 'CD4') ~ str_replace(source, 'CD4', 'CD4.'),
    str_detect(source, 'CD8') ~ str_replace(source, 'CD8', 'CD8.'),
    str_detect(source, 'CD14') ~ str_replace(source, 'CD14', 'CD14.'),
    str_detect(source, 'CD16') ~ str_replace(source, 'CD16', 'CD16.'),
    source == 'PreB' ~ 'Pre.B',
    source == 'LateEryth' ~ 'Late.Eryth',
    T ~ source
  ),
  target = case_when(
    target == 'HSCMPP' ~ 'HSC.MPP',
    str_detect(target, 'CD4') ~ str_replace(target, 'CD4', 'CD4.'),
    str_detect(target, 'CD8') ~ str_replace(target, 'CD8', 'CD8.'),
    str_detect(target, 'CD14') ~ str_replace(target, 'CD14', 'CD14.'),
    str_detect(target, 'CD16') ~ str_replace(target, 'CD16', 'CD16.'),
    target == 'PreB' ~ 'Pre.B',
    target == 'LateEryth' ~ 'Late.Eryth',
    T ~ target
  ))

# import omnipathdb to annotate interactions with curation effort
omnidb <- import_all_interactions() %>% dplyr::select(ligand = source_genesymbol, receptor = target_genesymbol, curation_effort, n_references, n_resources)
both_sig <- both_sig %>% left_join(omnidb, by = c('ligand', 'receptor'))

# filter to only interactions with curation effort >= 3
high_conf <- both_sig %>% filter(curation_effort >= 3) %>% 
  dplyr::select(source, target, ligand, receptor, aggregate_rank, sca.LRscore, tp, curation_effort)  

# define which timepoints each interaction was identified in 
high_conf <- high_conf %>% pivot_wider(c(source,target, ligand, receptor, curation_effort), names_from = 'tp',  values_from = c('aggregate_rank', 'sca.LRscore')) %>%
  mutate(timepoint = case_when(!is.na(aggregate_rank_Diagnosis) & !is.na(aggregate_rank_Healthy) ~ 'Both',
                               is.na(aggregate_rank_Diagnosis) ~ 'Healthy',
                               is.na(aggregate_rank_Healthy) ~ 'Diagnosis'))

# annotate each interaction with kegg pathways
prot <- unique(c(high_conf$ligand, high_conf$receptor))
msig <- import_omnipath_annotations(
  proteins = prot,
  resources = c('MSigDB'),
  wide = T
) 

kegg <- msig %>% filter(collection == 'kegg_pathways') %>% 
  dplyr::select(genesymbol, geneset) %>% 
  group_by(genesymbol) %>% 
  summarise(kegg_pathway=paste0(geneset,collapse = ", "))

annot_ints <- left_join(high_conf, kegg, by=c("ligand"="genesymbol")) %>% 
  left_join(kegg, by=c("receptor"="genesymbol"), suffix=c(".ligand", ".receptor"))

# add columns for number of cell types each interaction is identified in
annot_ints <- annot_ints %>% mutate(interacting_cell = ifelse(source == 'HSC.MPP', target, source),
                                    from_hsc = source == 'HSC.MPP',
                                    lr_pair = paste0(ligand, '|', receptor))
counts <- annot_ints %>% group_by(timepoint, lr_pair, from_hsc) %>% tally(name = 'number_of_celltypes')
annot_ints <- left_join(annot_ints, counts, by = c("timepoint", "lr_pair", "from_hsc"))
cts <- c()
for(i in 1:nrow(annot_ints)){
  cells <- annot_ints %>% filter(lr_pair == annot_ints$lr_pair[i], from_hsc == annot_ints$from_hsc[i], timepoint == annot_ints$timepoint[i]) %>% pull(interacting_cell) %>% unique() %>% paste0(collapse = ', ')
  cts <- c(cts, cells)
}
annot_ints$celltypes <- cts
annot_ints <- annot_ints %>% dplyr::select(interacting_cell, timepoint, source, target, ligand, receptor, aggregate_rank_Healthy, aggregate_rank_Diagnosis, sca.LRscore_Healthy, sca.LRscore_Diagnosis, curation_effort, kegg_pathway.ligand, kegg_pathway.receptor, number_of_celltypes, celltypes)

# save annotated interactions
saveRDS(annot_ints, 'annotated_interactions.Rds')

# for supplementary table
for(i in unique(annot_ints$interacting_cell)){
  sub <- annot_ints %>% filter(interacting_cell == i)
  sub <- sub %>% arrange(timepoint, aggregate_rank_Diagnosis, aggregate_rank_Healthy)
  fname <- paste0('int_csvs/', i, '.csv')
  write.table(sub, file = fname, quote = F, row.names = F, sep = '\t')
}
