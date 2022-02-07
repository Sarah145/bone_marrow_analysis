# Load libraries
library(tidyverse)
library(ggVennDiagram)
library(org.Hs.eg.db)
library(OmnipathR)
library(patchwork)
library(ggtext)
library(paletteer)

# Read in cell-cell interaction results and apply filters
healthy <- readRDS('~/Documents/Project/cc_interactions/healthy_sub1_liana_aggreg.Rds') %>% 
  filter(aggregate_rank <= 0.05, (source %in% c('HSCMPP') | target %in% c('HSCMPP'))) %>% 
  mutate(tp = 'Healthy')

d0 <- readRDS('~/Documents/Project/cc_interactions/d0_sub1_liana_aggreg.Rds') %>% 
  filter(aggregate_rank <= 0.05, (source %in% c('HSCMPP') | target %in% c('HSCMPP'))) %>% 
  mutate(tp = 'Diagnosis')
col_pals <- readRDS('colour_palettes.Rds')
cell_cols <- col_pals$celltype

# combine healthy & aml results
both_sig <- rbind(healthy, d0)

# annotate with curation effort using omnipath to get high confidence interactions
omnidb <- import_all_interactions() %>% dplyr::select(ligand = source_genesymbol, receptor = target_genesymbol, curation_effort, n_references, n_resources)
both_sig <- both_sig %>% left_join(omnidb, by = c('ligand', 'receptor'))
high_conf <- both_sig %>% filter(curation_effort >= 3)
high_conf <- high_conf %>% pivot_wider(c(source,target, ligand, receptor), names_from = 'tp',  values_from = c('aggregate_rank', 'sca.LRscore')) %>% 
  mutate(timepoint = case_when(!is.na(aggregate_rank_Diagnosis) & !is.na(aggregate_rank_Healthy) ~ 'Both',
                               is.na(aggregate_rank_Diagnosis) ~ 'Healthy',
                               is.na(aggregate_rank_Healthy) ~ 'Diagnosis'))

# Figure 4a
interactions <- list(Healthy = high_conf %>% filter(timepoint == 'Healthy' | timepoint == 'Both') %>% mutate(int = paste0(source, target, ligand, receptor)) %>% pull(int) %>% unique(),
                     `AML\n(diagnosis)` = high_conf %>% filter(timepoint == 'Diagnosis' | timepoint == 'Both') %>% mutate(int = paste0(source, target, ligand, receptor)) %>% pull(int) %>% unique())

ggVennDiagram(interactions, label_size = 12, set_size = 12, label_alpha = 1) +
  scale_color_manual(values = unname(c(col_pals$timepoint['Healthy'], col_pals$timepoint['AML (diagnosis)']))) +
  scale_fill_gradient(low = '#FFFFFF', high = '#FFFFFF') +
  theme(legend.position = 'none',
        text = element_text(size = 14))
ggsave('plots/venn.svg', width = 7, height = 5)

# Split common interactions into higher at healthy/AML
high_conf <- high_conf %>%  mutate(scaLogFC = log(sca.LRscore_Diagnosis/sca.LRscore_Healthy),
                                   tp = case_when(scaLogFC >= log(1.05) ~ 'Higher in AML',
                                                  scaLogFC <= log(0.95) ~ 'Higher in Healthy',
                                                  timepoint == 'Healthy' ~ 'Only in Healthy',
                                                  timepoint == 'Diagnosis' ~ 'Only in AML',
                                                  T ~ 'Unchanged')) 
# Annotate with KEGG and cellchat annotations
cellchat <- import_omnipath_annotations(
  proteins = unique(c(high_conf$ligand, high_conf$receptor)),
  resources = c('CellChatDB'),
  wide = TRUE
)

annot <- inner_join(
  high_conf, cellchat, # here we select the two dataframes
  by=c("ligand"="genesymbol")) %>% # the names of the columns to combine them on
  inner_join(
    cellchat, # the first dataframe is the result of the first operation here, second                              stays the same
    by=c("receptor"="genesymbol"), # the names of the columns to combine them on
    suffix=c(".ligand", ".receptor")) 

all_cellchat <- annot %>% filter(pathway.ligand == pathway.receptor) %>% 
  mutate(pathway = pathway.ligand) %>% dplyr::select(
    source, target, ligand, receptor, tp, pathway) %>%
  distinct() %>%
  mutate(other_cell = ifelse(source == 'HSCMPP', target, source),
         other_cell = case_when(
           other_cell == 'HSCMPP' ~ 'HSC.MPP',
           str_detect(other_cell, 'CD4') ~ str_replace(other_cell, 'CD4', 'CD4.'),
           str_detect(other_cell, 'CD8') ~ str_replace(other_cell, 'CD8', 'CD8.'),
           str_detect(other_cell, 'CD14') ~ str_replace(other_cell, 'CD14', 'CD14.'),
           str_detect(other_cell, 'CD16') ~ str_replace(other_cell, 'CD16', 'CD16.'),
           other_cell == 'PreB' ~ 'Pre.B',
           other_cell == 'LateEryth' ~ 'Late.Eryth',
           T ~ other_cell
         ),
         other_cell = factor(other_cell, levels = names(cell_cols))) %>%
  group_by(tp, other_cell, pathway) %>% tally() %>%
  mutate(type = 'CellChat')

kegg <- import_omnipath_annotations(
  proteins = unique(c(high_conf$ligand, high_conf$receptor)),
  resources = c('MSigDB'),
  wide = T
) %>% filter(collection == 'kegg_pathways')


annot <- inner_join(
  high_conf, kegg, # here we select the two dataframes
  by=c("ligand"="genesymbol")) %>% # the names of the columns to combine them on
  inner_join(
    kegg, # the first dataframe is the result of the first operation here, second                              stays the same
    by=c("receptor"="genesymbol"), # the names of the columns to combine them on
    suffix=c(".ligand", ".receptor")) 

all_kegg <- annot %>% filter(geneset.ligand == geneset.receptor) %>% 
  mutate(pathway = geneset.ligand) %>% dplyr::select(
    source, target, ligand, receptor, tp, pathway) %>%
  distinct() %>%
  mutate(pathway = str_sub(pathway, start = 6),
         pathway = str_replace_all(pathway, '_', ' '),
         other_cell = ifelse(source == 'HSCMPP', target, source),
         other_cell = case_when(
           other_cell == 'HSCMPP' ~ 'HSC.MPP',
           str_detect(other_cell, 'CD4') ~ str_replace(other_cell, 'CD4', 'CD4.'),
           str_detect(other_cell, 'CD8') ~ str_replace(other_cell, 'CD8', 'CD8.'),
           str_detect(other_cell, 'CD14') ~ str_replace(other_cell, 'CD14', 'CD14.'),
           str_detect(other_cell, 'CD16') ~ str_replace(other_cell, 'CD16', 'CD16.'),
           other_cell == 'PreB' ~ 'Pre.B',
           other_cell == 'LateEryth' ~ 'Late.Eryth',
           T ~ other_cell
         ),
         other_cell = factor(other_cell, levels = names(cell_cols))) %>%
  group_by(tp, other_cell, pathway) %>% tally() %>%
  mutate(type = 'KEGG')

all_pathways <- rbind(all_cellchat, all_kegg) %>% mutate(
  tp = factor(tp, levels = c('Only in Healthy', 'Higher in Healthy', 'Unchanged', 'Higher in AML', 'Only in AML')),
  type = factor(type, levels = c('KEGG', 'CellChat'))
) %>% group_by(pathway) %>% mutate(tot_n = sum(n)) %>% ungroup() %>% arrange(tot_n) %>% mutate(pathway = fct_inorder(pathway))

pathway_categories <- list(Hematopoiesis = c("NOTCH", "CD34", "KIT", "FLT3", "THBS", "HEMATOPOIETIC CELL LINEAGE", 
                                             "NOTCH SIGNALING PATHWAY", "TGF BETA SIGNALING PATHWAY"),
                           `Cytokine signalling` = c("TNF", "CXCL", "ICOS", "ADIPOCYTOKINE SIGNALING PATHWAY", "FLT3",         
                                                     "CHEMOKINE SIGNALING PATHWAY", "BAFF", "MK", "NATURAL KILLER CELL MEDIATED CYTOTOXICITY", 
                                                     "CYTOKINE CYTOKINE RECEPTOR INTERACTION"),
                           `Cancer-related` = c("SMALL CELL LUNG CANCER", "CHRONIC MYELOID LEUKEMIA", "COLORECTAL CANCER",       
                                                "PANCREATIC CANCER", "PATHWAYS IN CANCER"),
                           `Cell adhesion & homing` = c("CD34", "FN1", "CD99", "CLEC", "ADHERENS JUNCTION", "ICAM", "ALCAM", "CD6", "CD22",
                                                        "CD45", "REGULATION OF ACTIN CYTOSKELETON", "LEUKOCYTE TRANSENDOTHELIAL MIGRATION", 
                                                        "SELL", "FOCAL ADHESION", "ITGB2", "SELPLG", "ECM RECEPTOR INTERACTION", "ANTIGEN PROCESSING AND PRESENTATION", 
                                                        "CELL ADHESION MOLECULES CAMS"),
                           `Cellular stress & apoptosis` = c("TNF", "BAG", "APRIL", "ALZHEIMERS DISEASE", "APOPTOSIS", "APP", "MAPK SIGNALING PATHWAY")
)

int_types <- c()
for(i in 1:nrow(all_pathways)){
  int_types <- c(int_types, paste0(names(pathway_categories)[str_detect(pathway_categories, as.character(all_pathways$pathway[i]))], collapse = ', '))
}
all_pathways$int_type <- int_types

cat_labs <- c(`Cancer-related` = '<span style="color:#6EDCF9">\\*</span>',
              `Cell adhesion & homing` = '<span style="color:#FCE15E">\\*</span>',
              `Cellular stress & apoptosis` = '<span style="color:#E9CCD4">\\*</span>',
              `Cytokine signalling` = '<span style="color:#FD0409">\\*</span>',
              Hematopoiesis = '<span style="color:#471839">\\*</span>')
labs <- c()
for(i in 1:nrow(all_pathways)){
  l <- ifelse(all_pathways$int_type[i] == '', as.character(all_pathways$pathway[i]),
              paste0(paste0(cat_labs[str_split(as.character(all_pathways$int_type[i]), ', ')[[1]]], collapse = ''), all_pathways$pathway[i]))
  labs <- c(labs, l)
}
all_pathways$lab <- labs
ggplot(all_pathways %>% filter(other_cell != '**Total**', tp %in% c('Only in Healthy', 'Unchanged', 'Only in AML')) %>% 
         mutate(type = factor(type, levels = c('KEGG', 'CellChat')),
                int_type = factor(int_type, levels = names(cat_labs))) %>%
         arrange(int_type) %>%
         mutate(lab = fct_rev(fct_inorder(lab))), aes(x = n, y = lab, fill = other_cell)) +
  geom_col() +
  scale_fill_manual(values = cell_cols, name = 'Interacting cell') +
  scale_x_continuous(expand = c(0,0), name = '# interactions annotated') +
  labs(y = NULL) +
  facet_grid(type~tp, scales = 'free_y', space = 'free_y', switch = 'y') +
  theme_minimal(base_size = 20) +
  theme(axis.text = element_text(colour = 'black'),
        axis.text.y = element_markdown(),
        axis.ticks = element_line(colour = 'gray20'),
        strip.placement = 'outside',
        strip.text = element_text(face = 'bold', size = 20),
        legend.position = 'bottom',
        legend.box.margin = margin(r = -100),
        legend.text = element_text(size = 20))
ggsave('plots/functional_annotations.png', height = 15, width = 20, bg='white')

# Make legend
leg_df <- data.frame(`Category` = names(cat_labs), col = as.character(paletteer_d('fishualize::Parupeneus_insularis', 5)))
ggplot(leg_df, aes(x = Category, y = 5, col=  Category)) +
  geom_point(pch = '*', size = 15) +
  scale_color_paletteer_d('fishualize::Parupeneus_insularis') +
  guides(color = guide_legend(nrow = 3)) +
  theme_void(base_size = 30) +
  theme(legend.position = 'bottom',
        legend.key.height = unit(0, 'lines'))
ggsave('plots/pathway_cat_legend.png')
