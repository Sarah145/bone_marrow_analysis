# Load libraries
library(tidyverse)
library(ggtext)
library(ggrepel)
library(patchwork)

# Read in data 
df <- read.csv('~/Documents/Project/bone_marrow_analysis/data/integrated_umap3.csv', sep = '\t', row.names = 1) %>% 
  rownames_to_column() 

# Figure 1b
dataset_cols <- c(nuig_1 = '#4D1282', nuig_2 = '#4D1282', vanG = '#AD343E', petti = '#F2AF29', granja = '#B2EF9B', oetjen = '#6DD3CE')
ggplot(df, aes(x = umap1, y=umap2, col = dataset)) +
  geom_point(size = 0.1, show.legend = T, alpha = 0.6, shape='.') +
  scale_color_manual(values = dataset_cols, name = 'Dataset', breaks = c('nuig_1', 'granja', 'oetjen', 'petti', 'vanG'), labels = seq(1:5)) +
  guides(color = guide_legend(override.aes = list(shape = 19, size = 5, alpha = 1))) +
  theme_void(base_size = 20) +
  theme(legend.position = c(0.1, 0.8))
ggsave('plots/dataset_umap.png', width = 7, height = 7)

# Figure 1c
df <- df %>% mutate(timepoint = case_when(
  tp == 'Dia' ~ 'AML (diagnosis)',
  tp == 'Hea' ~ 'Healthy',
  tp == 'Rel' ~ 'AML (relapse)',
  tp == 'TRT' ~ 'AML (post-treatment)',
  tp == 'MRD' ~ 'AML (post-treatment)'
))
timepoint_cols <- c(Healthy = '#FF8965', `AML (diagnosis)` = '#9B7ED8', `AML (post-treatment)` = '#6C43C7', `AML (relapse)` = '#311C5E')
ggplot(df, aes(x = umap1, y=umap2, col = timepoint)) +
  geom_point(size = 0.1, show.legend = T, alpha = 0.6, shape='.') +
  scale_color_manual(values = timepoint_cols, name = 'Timepoint', labels = function(x) str_wrap(x, width = 10)) +
  guides(color = guide_legend(override.aes = list(shape = 19, size = 5, alpha = 1))) +
  theme_void(base_size = 20) +
  theme(legend.position = c(0.12, 0.8),
        legend.text = element_text(margin = margin(b = 2.5, t = 2.5)))
ggsave('plots/timepoint_umap.png', width = 7, height = 7)

# Figure 1d
cell_cols <- c(Plasma = "#333333", CLP = "#B2EF9B", Pre.B = "#3A7219", B = "#1B360C", 
               HSC.MPP = "#F7D760", MEP = "#B84F09", Late.Eryth = "#6D2F05", Megakaryocytes = "#F7A660", 
               GMP = "#DE9197", CD14.Mono = "#AD343E", CD16.Mono = "#BD0000", cDC = "#FF0378", 
               pDC = "#FC5E03", CD4.N = "#015073", CD4.M = "#96CBFE", CD8.N = "#16B8B2", 
               CD8.M = "#946CD5", NK = "#4D1282", Mesenchymal = "#702127")

df$ct <- factor(df$celltype_final, levels = names(cell_cols))
ggplot(df, aes(x = umap1, y=umap2, col = ct)) +
  geom_point(size = 0.1, show.legend = F, alpha = 0.6, shape='.') +
  geom_label_repel(data = df %>%
                     group_by(ct) %>%
                     summarise(x = median(umap1),
                               y = median(umap2)),
                   aes(x = x,
                       y = y,
                       label = levels(df$ct)), 
                   show.legend = F,
                   label.size = NA,
                   label.padding = unit(0.1, "lines"),
                   fill = alpha(c("white"),0.8),
                   size = 6, max.overlaps = 50,
                   nudge_y = 0.5, segment.colour = NA) +
  scale_color_manual(values = cell_cols, name = NULL) +
  guides(color = guide_legend(override.aes = list(shape = 19, size = 3, alpha = 1))) +
  theme_void(base_size = 14) 
ggsave('plots/celltype_umap.svg', width = 7, height = 7) # save as svg so I can move labels in inkscape

# Figure 1e
# Read in marker gene expression data
gene_df <- read.csv('~/Documents/Project/bone_marrow_analysis/data/marker_gene_exp1.csv.gz', sep = '\t', row.names = 1) %>% 
  rownames_to_column() 

df <- left_join(df, gene_df, by = 'rowname')  

# get scaled average expression of each gene in each cell type and
# percent of each cell type expressing each gene
gen_df <- df %>% select(14:ncol(df)) %>% 
  pivot_longer(!ct, names_to = 'gene') %>% 
  group_by(ct, gene) %>% 
  mutate(per_exp = sum(value > 0)/length(ct), 
         avg_exp = mean(value)) %>% 
  select(-value) %>% 
  distinct() %>% ungroup() %>% group_by(gene) %>% 
  mutate(scaled_exp = avg_exp/max(avg_exp))

gene_order <- c('SDC1', 'MME', 'CD19', 'CD79A', 'MS4A1', 'KIT', 'CD34', 'CLC',
                'GYPA', 'PPBP', 'ITGA2B',  'CD33', 'MPO', 'CTSG', 'AZU1', 'LYZ', 'S100A8', 'S100A9', 'CD14', 'ITGAM', 'FCGR3A',
                'FCER1G', 'FCER1A', 'CLEC10A', 'IL3RA', 'IRF8',
                'CD4', 'IL7R', 'CD3E', 'CCR7', 'CD8A', 'CCL5', 'GNLY', 'KLRB1', 'NKG7',
                'CXCL12', 'VCAM1')

grad_cols <- colorRampPalette(c('#F7D340', '#E55C30', '#87216B', '#190C3E'))
y_labs <- paste0("<span style='color:", cell_cols, "'>", names(cell_cols), "</span>")

ggplot(gen_df %>% filter(gene %in% gene_order), 
       aes(x = gene, y = ct, size = per_exp, col = scaled_exp)) +
  geom_point() +
  scale_y_discrete(limits = rev(names(cell_cols)), labels = rev(y_labs)) +
  scale_x_discrete(limits = gene_order) +
  scale_colour_gradientn(colours = rev(grad_cols(100))) +
  scale_size(labels = scales::percent, range = c(0.5, 4)) +
  labs(x = NULL, y = NULL, colour = 'Avg. expression\n(scaled 0-1)', size = '% cells\nexpressing\ngene') +
  guides(colour = guide_colourbar(order = 1, barlength = unit(6, 'lines')), size = guide_legend(order = 2)) +
  theme_linedraw(base_size = 14) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.text.y = element_markdown(face = 'bold'),
        legend.title = element_text(size = 12, vjust = 0.5),
        legend.box.margin = margin(t = 45)
        )
ggsave('plots/marker_genes.png', width = 10, height = 4)

# Save colour palettes to use for other plots
col_pals <- list(dataset = dataset_cols, timepoint = timepoint_cols, celltype = cell_cols)
saveRDS(col_pals, 'colour_palettes.Rds')
