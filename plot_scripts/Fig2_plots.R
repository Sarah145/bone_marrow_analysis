# Load libraries
library(tidyverse)
library(ggtext)
library(patchwork)
library(cowplot)
library(gridExtra)

# Figure 2a
# Read in data
df <- read.csv('~/Documents/Project/bone_marrow_analysis/data/integrated_umap3.csv', sep = '\t', row.names = 1) %>%
  filter(sample %ni% c('CD34_D2T1', 'CD34_D3T1', 'Normal_sorted_170531', 'Normal_sorted_170607', 'BM5-34p', 'BM5-34p38n')) %>% # remove sorted samples
  mutate(timepoint = case_when(
    tp == 'Dia' ~ 'AML (diagnosis)',
    tp == 'Hea' ~ 'Healthy',
    tp == 'Rel' ~ 'AML (relapse)',
    tp == 'TRT' ~ 'AML (post-treatment)',
    tp == 'MRD' ~ 'AML (post-treatment)'
  ))

df$type <- case_when(
  df$celltype_final %in% c('CD4.M', 'CD4.N', 'CD8.M', 'CD8.N', 'CLP', 'Pre.B', 'B', 'Plasma', 'NK') ~ 'Lymphoid',
  df$celltype_final %in% c('HSC.MPP', 'GMP', 'cDC', 'pDC', 'CD14.Mono', 'CD16.Mono') ~ 'Myeloid',
  df$celltype_final %in% c('MEP', 'Late.Eryth') ~ 'Erythroid',
  df$celltype_final %in% c('Mesenchymal', 'Megakaryocytes') ~ 'Other')

ggplot(df, aes(x = timepoint, fill = type)) +
  geom_bar(position = 'fill', show.legend = T) +
  scale_x_discrete(limits = c('Healthy', 'AML (diagnosis)', 'AML (post-treatment)', 'AML (relapse)'), labels = function(x) str_wrap(x, width = 10)) +
  scale_y_continuous(labels = scales::percent, expand = c(0,0)) +
  scale_fill_manual(values =c(Myeloid = '#00ccff', Lymphoid = '#ff8a44', Erythroid = '#ee4952', Other = '#90005d'), name = NULL) +
  theme_classic(base_size = 23) +
  theme(axis.text = element_text(colour = 'black'),
        axis.title = element_blank(),
        panel.background = element_rect(colour = 'transparent', fill = 'transparent'),
        plot.background = element_rect(colour = 'transparent', fill = 'transparent'))
ggsave('plots/type_timepoint_comp.png', width = 10, height = 6, bg = 'transparent')

cell_props <- readRDS('../propeller_results.Rds')
cell_props <-  cell_props %>%
  mutate(timepoint = case_when(
    tp == 'AML.Diagnosis' ~ 'AML (diagnosis)',
    tp == 'Healthy' ~ 'Healthy',
    tp == 'AML.Relapse' ~ 'AML (relapse)',
    tp == 'AML.Posttreatment' ~ 'AML (post-treatment)'
  ),
  TP1 = case_when(
    TP1 == 'AML.Diagnosis' ~ 'AML (diagnosis)',
    TP1 == 'Healthy' ~ 'Healthy',
    TP1 == 'AML.Relapse' ~ 'AML (relapse)',
    TP1 == 'AML.Posttreatment' ~ 'AML (post-treatment)'
  ),
  TP2 = case_when(
    TP2 == 'AML.Diagnosis' ~ 'AML (diagnosis)',
    TP2 == 'Healthy' ~ 'Healthy',
    TP2 == 'AML.Relapse' ~ 'AML (relapse)',
    TP2 == 'AML.Posttreatment' ~ 'AML (post-treatment)'
  ))

col_pals <- readRDS('colour_palettes.Rds')
cell_cols <- col_pals$celltype

cell_props$celltype <- factor(cell_props$BaselineProp.clusters, levels = names(cell_cols))
cell_props$timepoint <- factor(cell_props$timepoint, levels = c('Healthy', 'AML (diagnosis)', 'AML (post-treatment)', 'AML (relapse)'))
cell_props$TP1 <- factor(cell_props$TP1, levels = c('Healthy', 'AML (diagnosis)', 'AML (post-treatment)', 'AML (relapse)'))
cell_props$TP2 <- factor(cell_props$TP2, levels = c('Healthy', 'AML (diagnosis)', 'AML (post-treatment)', 'AML (relapse)'))




df$celltype <- factor(df$celltype_final, levels = names(cell_cols))
plot_df <- df %>% group_by(sample, timepoint, celltype) %>% tally() %>% 
  ungroup() %>% group_by(sample) %>% mutate(per = n/sum(n))

plot_list <- list()
for(i in levels(cell_props$celltype)){
  sub <- cell_props %>% filter(celltype == i)
  plot_df_sub <- plot_df %>% filter(celltype == i) 
  bar_plot <- ggplot(sub %>% select(timepoint, celltype, prop) %>% mutate(prop = round(prop, 5)) %>% distinct(),
                     aes(x = timepoint, y = prop, fill = celltype)) +
    geom_col(show.legend = F) +
    geom_jitter(data = plot_df_sub, aes(x = timepoint, y = per, col = NULL, fill = NULL), size = 2, alpha = 0.5, position = position_jitter(width = 0.1), show.legend = F) +
    scale_fill_manual(values = cell_cols) +
    scale_x_discrete(labels = function(x) str_wrap(x, width = 10)) +
    scale_y_continuous(expand = c(0,0), limits = c(0, max(plot_df_sub$per) * 1.1)) +
    facet_wrap(~celltype, scales = 'free') +
    labs(y = 'Proportion', x = NULL) +
    theme_classic(base_size = 14) +
    theme(axis.text = element_text(colour = 'black'),
          axis.text.x = element_text(margin = margin(t = 10, b = -5), vjust = 0.5),
          axis.ticks = element_line(colour = 'gray20'),
          strip.text = element_text(colour = 'black', size = 13, face = 'bold'))
  
  d <- data.frame()
  for(j in levels(sub$TP1)){
    for(k in levels(sub$TP1)){
      sub_sub <- sub %>% filter((TP1 == j & TP2 == k) | (TP2 == j & TP1 == k))
      enriched <- F
      if(nrow(sub_sub) > 0){
        if(j == unique(sub_sub$TP1) & k == unique(sub_sub$TP2) & unique(sub_sub$TP1_greater)){enriched <- T}
        if(k == unique(sub_sub$TP1) & j == unique(sub_sub$TP2) & unique(sub_sub$TP1_greater) == F){enriched <- T}}
      dd <- data.frame(TP1 = j, TP2 = k, enr = enriched, fdr = ifelse(nrow(sub_sub) > 0, unique(sub_sub$p_adj), NA), celltype = i)
      d <- rbind(d, dd)
    }
  }
  
  tile_plot <- ggplot(d,
                      aes(x = TP1, y = TP2)) +
    geom_tile(fill = 'white', col = 'gray20', size = 0.5) +
    geom_point(data = d %>% filter(fdr <= 0.05, enr == T),aes(size = -log10(fdr)), show.legend = F, col = cell_cols[i]) +
    scale_x_discrete(limits = levels(cell_props$TP1), position = 'top', expand = c(0,0)) +
    scale_y_discrete(limits = rev(levels(cell_props$TP2)), expand = c(0,0), labels = function(x) str_wrap(x, width = 10)) +
    scale_size(limits = c(0,20), breaks = seq(5,20,5), range = c(0.5, 6)) +
    facet_wrap(~celltype, scales = 'free') +
    labs(y = 'Compared to') +
    theme_linedraw(base_size = 14) +
    theme(axis.text = element_text(colour = 'black'),
          axis.title.x = element_blank(),
          axis.text.y = element_text(hjust = 0.5),
          axis.ticks = element_line(colour = 'gray20'),
          strip.text = element_blank(),
          axis.text.x.top = element_blank())
  comb_plot <- bar_plot/tile_plot
  plot_list[[i]] <- eval(substitute(comb_plot))
}

g <- plot_grid(plotlist = plot_list, align = 'hv', nrow = 4)
ggdraw(g) + 
  theme(plot.background = element_rect(fill="white", color = NA))

ggsave('plots/propeller_results.png', width = 24, height = 20)


dummy_df <- data.frame(enriched = c(T, F), fdr = c(1, 20))

ggplot(dummy_df %>% mutate(Enriched = factor(enriched, levels = c(T, F))), aes(x = 1, y = 2)) +
  geom_tile(aes(col = Enriched), fill = 'white', size = 0.5) +
  geom_point(aes(col = Enriched, alpha = Enriched, size = fdr)) +
  scale_alpha_manual(values = c(1,0)) +
  scale_colour_manual(values = c('black', 'black')) +
  scale_size(limits = c(0,20), breaks = seq(5,20,5), range = c(0.5, 6), name = '-log10(FDR)') +
  theme_void(base_size = 20) +
  theme(legend.direction = 'vertical',
        legend.position = 'bottom',
        legend.spacing = unit(2, 'lines'),
        plot.margin = margin(b = 10))
ggsave('~/Documents/Project/my_first_paper_<3/fig2_legend.png', width = 3, height = 2)


# ARC talk plot
ggplot(cell_props %>% select(timepoint, celltype, prop) %>% mutate(prop = round(prop, 5)) %>% distinct(),
       aes(x = timepoint, y = prop, fill = celltype)) +
  geom_col(show.legend = F) +
  geom_jitter(data = df1, aes(x = timepoint, y = per, col = NULL, fill = NULL), size = 0.2, alpha = 0.5, position = position_jitter(width = 0.1), show.legend = F) +
  scale_fill_manual(values = cell_cols) +
  scale_x_discrete(labels = function(x) str_wrap(x, width = 10)) +
  scale_y_continuous(expand = c(0,0)) +
  facet_wrap(~celltype, scales = 'free') +
  labs(y = 'Proportion', x = NULL) +
  theme_classic(base_size = 18) +
  theme(axis.text = element_text(colour = 'black', size = 16),
        axis.text.x = element_text(margin = margin(b = -5), vjust = 0.5, size = 16, lineheight = 0.5),
        axis.ticks = element_line(colour = 'gray20', size = 0.2),
        axis.line = element_line(size = 0.4),
        strip.background = element_rect(size = 0.5, fill = 'transparent'),
        strip.text = element_text(colour = 'black', size = 18, face = 'bold')) +
  coord_cartesian(clip= 'off')
ggsave('~/Documents/Project/pics/propeller_results_bar.png', width = 10, height = 5)
