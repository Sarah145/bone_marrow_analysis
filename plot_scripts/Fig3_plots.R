# Load libraries
library(tidyverse)
library(ggtext)
library(paletteer)

# Read in DEG and GO results
degs <- readRDS('../data/d0_healthy_degs_libra.Rds')
go_df <- readRDS('../data/healthy_d0_enriched_gos_1FC.Rds')

# Figure 3a
col_pals <- readRDS('colour_palettes.Rds')
cell_cols <- col_pals$celltype
lab <- paste0("<span style='color:", cell_cols, "'>", names(cell_cols), "</span>")
names(lab) <- names(cell_cols)
degs$cell_type <- factor(degs$cell_type, levels = names(cell_cols))
cols <- paletteer_d('unikn::pal_signal', 3)
names(cols) <- c('Downregulated', 'NS', 'Upregulated')
degs <- degs %>% mutate(type = case_when(
  avg_logFC >= 1 & p_val_adj <= 0.05 ~ 'Upregulated',
  avg_logFC <= -1 & p_val_adj <= 0.05 ~ 'Downregulated',
  T ~ 'NS'),
  cell_type = factor(cell_type, levels = names(cell_cols)))
label_df <- degs %>% group_by(cell_type, type) %>% tally() %>% filter(type != 'NS')

ggplot(degs, aes(x = avg_logFC, y = -log10(p_val_adj), col = type, shape = sig)) +
  geom_point(show.legend = T, size = 0.5) +
  geom_text(data = label_df %>% filter(type == 'Upregulated'), aes(x = 5, y =60, label = paste0('n=',n), col = type, shape = NULL)) +
  geom_text(data = label_df %>% filter(type == 'Downregulated'), aes(x = -5, y =60, label = paste0('n=',n), col = type, shape = NULL)) +
  scale_shape_manual(values = c(20,19)) +
  scale_color_manual(values = cols, name = NULL, breaks = c('Upregulated', 'Downregulated', 'NS')) +
  labs(x = 'Avg. logFC', y = '-log10(adj. p val)') +
  guides(colour = guide_legend(override.aes = list(size = 3)), shape = 'none') +
  facet_wrap(~cell_type, nrow = 2) +
  theme_minimal(base_size = 14) +
  theme(strip.background = element_rect(colour = 'black'),
        strip.text = element_markdown(face = 'bold', colour = 'black'),
        axis.text = element_text(colour = 'black'),
        legend.position = c(0.95, 0.17),
        legend.background = element_rect(colour = 'black'))  
ggsave('plots/volcanos.png', width = 18, height = 4.2)


# Figure 3b
top_degs <- degs %>% filter(sig == T) %>% mutate(cell_type = factor(cell_type, levels = names(cell_cols))) %>% group_by(cell_type) %>% slice_min(order_by = p_val_adj, n = 10) 
lab <- paste0("<span style='color:", cell_cols, "'>■</span>")
ggplot(top_degs %>% mutate(gene = fct_inorder(gene)), aes(y = cell_type, x = gene, size = -log10(p_val_adj), col = avg_logFC)) +
  geom_point() +
  scale_y_discrete(limits = rev(names(cell_cols)), labels = rev(lab)) +
  scale_x_discrete(position = 'top') +
  scale_size(range = c(0.5, 4)) +
  scale_colour_paletteer_c(palette = 'ggthemes::Classic Red-Blue') +
  labs(x = NULL, y = NULL, col = 'Avg. logFC', size = '-log10(adj. p val)') +
  guides(colour = guide_colourbar(order = 1)) +
  theme_minimal(base_size = 14) +
  theme(axis.text.x.top = element_text(angle = 90, hjust = 0, vjust = 0.5, colour = 'black', size = 10.5),
        axis.text.y = element_markdown(face = 'bold', size = 15, angle = 90, hjust = 0.5, vjust = 0),
        axis.ticks = element_line(colour = 'gray20'),
        legend.position = 'bottom') +
  coord_cartesian(clip = 'off')
ggsave('plots/top_deg_dotplot.png', width = 14, height = 5)


# Figure 3c
plot_go <- go_df %>% group_by(cell_type, timepoint) %>% filter(!is.na(qvalue), qvalue <= 0.05) %>% slice_min(order_by = qvalue, n= 5, with_ties = F)
grad_cols <- colorRampPalette(c('#F7D340', '#E55C30', '#87216B', '#190C3E'))

ggplot(plot_go %>% ungroup() %>% mutate(Description = fct_rev(fct_inorder(Description))), aes(x = cell_type, y = Description, size = -log10(qvalue), col = Count)) +
  geom_point() +
  scale_size(range = c(0.5, 4)) +
  scale_x_discrete(limits = names(cell_cols), labels = lab, position = 'top') +
  scale_colour_gradientn(colours = grad_cols(100)) +
  facet_wrap(~timepoint, scales = 'free_y') +
  theme_minimal(base_size = 14) +
  theme(axis.text.x.top = element_markdown(face= 'bold', size = 12),
        axis.text.y = element_text(colour = 'black'),
        axis.title.x = element_blank(),
        strip.placement = 'outside',
        strip.text = element_blank(),
        axis.ticks = element_line(colour = 'gray20'),
        legend.position = 'bottom')
ggsave('plots/top_gos.png', width = 14, height = 6)
