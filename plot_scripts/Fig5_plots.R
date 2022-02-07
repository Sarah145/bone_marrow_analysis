# Load libraries
library(tidyverse)
library(ggtext)
library(patchwork)

int_df <- readRDS('~/Documents/Project/my_first_paper_<3/annotated_interactions.Rds')
col_pals <- readRDS('colour_palettes.Rds')
cell_cols <- col_pals$celltype

# Figure 5a 
cyto <- int_df %>% filter(str_detect(kegg_pathway.ligand, 'CYTOKINE'), str_detect(kegg_pathway.receptor, 'CYTOKINE')) %>% 
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
  ),
  source = factor(source, levels = names(cell_cols)),
  target = factor(target, levels = names(cell_cols)),
  from_hsc = factor(source == 'HSC.MPP', levels = c(T,F), labels = c('From HSC.MPP', 'To HSC.MPP')),
  interacting_cell = case_when(from_hsc == 'From HSC.MPP' ~ target,
                               from_hsc == 'To HSC.MPP' ~ source),
  interacting_cell = factor(interacting_cell, levels = names(cell_cols))
  ) %>%
  arrange(from_hsc, source, target) %>% 
  pivot_longer(c(agg_rank_Healthy, agg_rank_Diagnosis), names_to = 'tp', values_to = 'score') %>%
  mutate(tp = factor(ifelse(tp == 'agg_rank_Healthy', 'Healthy', 'AML (diagnosis)'), levels = c('Healthy', 'AML (diagnosis)')),
         lr_pair = paste0(ligand, '|', receptor)) %>%
  arrange(lr_pair) %>%
  mutate(lr_pair = fct_rev(fct_inorder(lr_pair)))


cols <- colorRampPalette(c('#5CC8FF', '#E0479E'))

labs <- paste0("<span style='color:", cell_cols, "'>", names(cell_cols), "</span>")

ggplot(cyto, aes(x=interacting_cell, y = lr_pair,
                 col = -log10(score), size = -log10(score))) +
  geom_point() +
  scale_colour_gradientn(colours = cols(100), breaks = seq(2,4,1)) +
  scale_size(breaks = seq(2,4,1), range = c(0.5, 5)) +
  scale_x_discrete(limits = names(cell_cols), labels = labs) +
  labs(y = 'Ligand|Receptor', x = 'Interacting cell', size = 'Interaction\nscore', title = 'CYTOKINE CYTOKINE RECEPTOR INTERACTIONS') +
  facet_grid(from_hsc ~ tp, scales = 'free_y', switch = 'y', space = 'free_y') +
  guides(size = guide_legend(override.aes = list(colour = c("#5CC8FF", "#9E87CE", "#E0479E"))), colour = 'none') +
  theme_minimal(base_size = 14) +
  theme(strip.placement = 'outside',
        strip.text = element_text(colour = 'black'),
        strip.background = element_rect(colour = 'black'),
        axis.text = element_text(colour='black'),
        axis.text.x = element_markdown(face='bold',angle = 90, hjust = 1, vjust = 0.5),
        axis.ticks = element_line(colour='gray20'),
        plot.title = element_text(hjust = 0.5)
  ) +
  coord_cartesian(clip = 'off')
ggsave('plots/cytokine_interactions.png', width = 12, height = 7)


# Figure 5b
cellcell <- int_df %>% filter(str_detect(kegg_pathway.ligand, 'CELL_ADHESION'), str_detect(kegg_pathway.receptor, 'CELL_ADHESION')) %>% 
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
  ),
  source = factor(source, levels = names(cell_cols)),
  target = factor(target, levels = names(cell_cols)),
  from_hsc = factor(source == 'HSC.MPP', levels = c(T,F), labels = c('From HSC.MPP', 'To HSC.MPP')),
  interacting_cell = case_when(from_hsc == 'From HSC.MPP' ~ target,
                               from_hsc == 'To HSC.MPP' ~ source),
  interacting_cell = factor(interacting_cell, levels = names(cell_cols))
  ) %>%
  arrange(from_hsc, source, target) %>% 
  pivot_longer(c(agg_rank_Healthy, agg_rank_Diagnosis), names_to = 'tp', values_to = 'score') %>%
  mutate(tp = factor(ifelse(tp == 'agg_rank_Healthy', 'Healthy', 'AML (diagnosis)'), levels = c('Healthy', 'AML (diagnosis)')),
         lr_pair = paste0(ligand, '|', receptor)) %>%
  arrange(lr_pair) %>%
  mutate(lr_pair = fct_rev(fct_inorder(lr_pair)))

ggplot(cellcell, aes(x=interacting_cell, y = lr_pair,
                     col = -log10(score), size = -log10(score))) +
  geom_point() +
  scale_colour_gradientn(colours = cols(100), breaks = seq(2,5,1)) +
  scale_size(breaks = seq(2,5,1), range = c(0.5, 5)) +
  scale_x_discrete(limits = names(cell_cols), labels = labs) +
  labs(y = 'Ligand|Receptor', x = 'Interacting cell', size = 'Interaction\nscore', title = 'CELL ADHESION MOLECULES CAMS') +
  facet_grid(from_hsc ~ tp, scales = 'free_y', switch = 'y', space = 'free_y') +
  guides(size = guide_legend(override.aes = list(colour = c("#5CC8FF", "#889DDE", "#B372BE", "#E0479E"))), colour = 'none') +
  theme_minimal(base_size = 14) +
  theme(strip.placement = 'outside',
        strip.text = element_text(colour = 'black'),
        strip.background = element_rect(colour = 'black'),
        axis.text = element_text(colour='black'),
        axis.text.x = element_markdown(face='bold',angle = 90, hjust = 1, vjust = 0.5),
        axis.ticks = element_line(colour='gray20'),
        plot.title = element_text(hjust = 0.5)
  ) +
  coord_cartesian(clip = 'off')
ggsave('plots/cell_adhesion_interactions.png', width = 12, height = 7.5)
