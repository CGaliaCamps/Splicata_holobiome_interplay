rm(list = ls())

library(data.table)
library(arrow)
library(ggplot2)
library(ggsignif)
library(ggpubr)
library(stringr)
library(tidyr)
library(dplyr)
library(ggtext)
library(ggrepel)
library(tidyverse)

setwd("E:/Splicata/16S-RNA-DNA/Micro_function")

##########################################################################################
corr.gene.dbb = read.csv("KEGG_CAZyme_correlated_genome_dbb.csv", header = T)
length(unique(corr.gene.dbb$Genome)) # 19 genomes, correct
head(corr.gene.dbb)
indic.gene.dbb = read.csv("KEGG_CAZyme_indicator_genome_dbb.csv", header = T)
length(unique(indic.gene.dbb$Genome)) # 32 genomes, correct
head(indic.gene.dbb)

# 1 KEGG subsets
kegg.target.genes = read.csv("KEGG_symbiont_pathogen_expanded.csv", header = T)
length(unique(kegg.target.genes$KEGG)) # 150 KEGGs
head(kegg.target.genes)

corr.gene.dbb.kegg = subset(corr.gene.dbb, Database == "KOfam")
names(corr.gene.dbb.kegg)[names(corr.gene.dbb.kegg) == "Accession"] = "KEGG"
corr.gene.dbb.kegg.target = corr.gene.dbb.kegg[corr.gene.dbb.kegg$KEGG %in% kegg.target.genes$KEGG,]
corr.gene.dbb.kegg.target = merge(corr.gene.dbb.kegg.target, kegg.target.genes, by = "KEGG")
head(corr.gene.dbb.kegg.target)
length(unique(corr.gene.dbb.kegg.target$KEGG)) #76/150 KEGGs in our genomes
length(unique(corr.gene.dbb.kegg.target$Genome)) # 19/19 genomes contain some of these KEGGs

indic.gene.dbb.kegg = subset(indic.gene.dbb, Database == "KOfam")
names(indic.gene.dbb.kegg)[names(indic.gene.dbb.kegg) == "Accession"] = "KEGG"
indic.gene.dbb.kegg.target = indic.gene.dbb.kegg[indic.gene.dbb.kegg$KEGG %in% kegg.target.genes$KEGG,]
indic.gene.dbb.kegg.target = merge(indic.gene.dbb.kegg.target, kegg.target.genes, by = "KEGG")
head(indic.gene.dbb.kegg.target)
length(unique(indic.gene.dbb.kegg.target$KEGG)) #84/150 KEGGs in our genomes
length(unique(indic.gene.dbb.kegg.target$Genome)) # 32/32 genomes contain some of these KEGGs

# Plot frequencies of pathways - KEGG + Correlated
# Summarise: count KEGGs per Species, Category and get Correlation and Lifestyle
df_heatmap <- corr.gene.dbb.kegg.target %>%
  group_by(Species, Category, Correlation, Lifestyle) %>%
  summarise(n = n(), .groups = "drop")

# Order Species by Correlation for grouping on the axis
species_order <- corr.gene.dbb.kegg.target %>%
  select(Species, Correlation) %>%
  distinct() %>%
  arrange(Correlation, Species) %>%
  pull(Species)

# Order Categories by Lifestyle for grouping on the axis
category_order <- corr.gene.dbb.kegg.target %>%
  select(Category, Lifestyle) %>%
  distinct() %>%
  group_by(Lifestyle) %>%
  arrange(Lifestyle, desc(Category)) %>%  # reverse alphabetical *within* Lifestyle
  ungroup() %>%
  pull(Category)

df_heatmap <- df_heatmap %>%
  mutate(
    Species  = factor(Species, levels = unique(species_order)),
    Category = factor(Category, levels = unique(category_order))
  )

# Colour palette for Lifestyle (for Category axis labels)

lifestyle_colors <- c(
  "Symbiont"      = "#D85A30",
  "Pathogen"      = "#534AB7",
  "Both"          = "#2F7E36")

# Build a lookup for colouring Category axis text
category_lifestyle <- corr.gene.dbb.kegg.target %>%
  select(Category, Lifestyle) %>%
  distinct() %>%
  arrange(Lifestyle, Category)

category_label_df <- corr.gene.dbb.kegg.target %>%
  select(Category, Lifestyle) %>%
  distinct() %>%
  mutate(
    label_colored = paste0(
      "<span style='color:",
      lifestyle_colors[Lifestyle],
      "'>",
      Category,
      "</span>"
    )
  )

axis_colors <- lifestyle_colors[category_lifestyle$Lifestyle]
names(axis_colors) <- category_lifestyle$Category
axis_colors_ordered <- axis_colors[levels(df_heatmap$Category)]

df_heatmap.chronly = df_heatmap[df_heatmap$Correlation %in% c("pop", "chr11_onlychr", "chr2_onlychr"),]
df_heatmap.all = df_heatmap[df_heatmap$Correlation %in% c("pop", "chr11", "chr2", "chr4"),]

df_heatmap.chronly <- df_heatmap.chronly %>%
  mutate(Correlation = recode(
    Correlation,
    pop = "pop",
    chr11_onlychr = "chr11 only",
    chr2_onlychr  = "chr2 only"
  ))

# Heatmap - facets
f_heatmap <- ggplot(df_heatmap.chronly, aes(x = Species, y = Category, fill = n)) +
  geom_tile(colour = "white", linewidth = 0.3) +
  scale_fill_gradient(
    low = "#EAF3DE", high = "#173404",
    name = "KEGG count", na.value = "grey95"
  ) +
  scale_y_discrete(
    drop = FALSE,
    labels = setNames(
      category_label_df$label_colored,
      category_label_df$Category
    )
  ) +
  facet_grid(
    . ~ Correlation,
    scales = "free_x",
    space = "free_x",
    labeller = labeller(
      Correlation = function(x) stringr::str_wrap(x, 12)
    )
  ) +
  labs(
    x = NULL, y = NULL,
    title = "KEGG repertoire per correlated genome"
  ) +
  theme_bw(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
    axis.text.y = element_markdown(size = 12),
    strip.text = element_text(face = "bold", size = 10),
    strip.background = element_rect(fill = "grey92"),
    legend.position = "right",
    legend.key.size = unit(0.5, "cm"),
    panel.grid = element_blank(),
    plot.title = element_text(size = 14, face = "bold")
  )
f_heatmap

tiff("KEGG_heatmap_correlated_genomes_chromosome_only_2026.tiff", units="cm", width=30, height=16, res = 300, compression = 'lzw')
f_heatmap
dev.off()


# Heatmap with sorting by number of KEGG categories per genome
# Order Species by number of Categories present
species_order_by_categories <- df_heatmap %>%
  distinct(Species, Category) %>%     # ensure 1 row per Species-Category
  count(Species, name = "n_categories") %>%
  arrange(desc(n_categories)) %>%
  pull(Species)

# Apply ordering
df_heatmap <- df_heatmap %>%
  mutate(
    Species = factor(Species, levels = species_order_by_categories)
  )

n_heatmap <- ggplot(df_heatmap, aes(x = Species, y = Category, fill = n)) +
  geom_tile(colour = "white", linewidth = 0.3) +
  scale_fill_gradient(
    low = "#EAF3DE", high = "#173404",
    name = "KEGG count", na.value = "grey95"
  ) +
  scale_y_discrete(
    drop = FALSE,
    labels = setNames(
      category_label_df$label_colored,
      category_label_df$Category
    )
  ) +
  labs(
    x = NULL, y = NULL,
    title = "KEGG repertoire per correlated genome"
  ) +
  theme_bw(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14, face = "italic"),
    axis.text.y = element_markdown(size = 12),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 10),
    legend.position = "right",
    legend.key.size = unit(0.5, "cm"),
    panel.grid = element_blank(),
    plot.title = element_text(size = 14, face = "bold")
  )
n_heatmap

tiff("KEGG_heatmap_correlated_genomes_nofacet_2026.tiff", units="cm", width=24, height=16, res = 300, compression = 'lzw')
n_heatmap
dev.off()


# Plot as ordination
df_heatmap = as.data.frame(df_heatmap)

kegg_mat <- df_heatmap %>%
  as.data.frame() %>%
  group_by(Species, Category) %>%
  summarise(n = sum(n), .groups = "drop") %>%
  tidyr::pivot_wider(
    names_from  = Category,
    values_from = n,
    values_fill = list(n = 0)
  ) %>%
  tibble::column_to_rownames("Species")

pca <- prcomp(kegg_mat, center = TRUE, scale. = TRUE)
summary(pca)$importance[2, 1:2]
scores <- as.data.frame(pca$x[, 1:2])
scores$Species <- rownames(scores)
loadings <- as.data.frame(pca$rotation[, 1:2])
loadings$Category <- rownames(loadings)

arrow_scale <- 2.5
loadings_plot <- loadings %>%
  mutate(
    PC1 = PC1 * arrow_scale,
    PC2 = PC2 * arrow_scale
  )
loadings_plot <- loadings_plot %>%
  filter(abs(PC1) + abs(PC2) > 0.9)

lifestyle_colors <- c(
  Symbiont  = "#D85A30",
  Pathogen  = "#534AB7",
  Both      = "#2F7E36"
)

category_lifestyle <- df_heatmap %>%
  select(Category, Lifestyle) %>%
  distinct()

loadings_plot <- loadings_plot %>%
  left_join(category_lifestyle, by = "Category")

scores_labels <- scores %>%
  filter(PC1 < 0)

pca.corr.kegg <- ggplot() +
  
  # Genome points
  geom_point(
    data = scores,
    aes(PC1, PC2),
    size = 3.5,
    color = "black"
  ) +
  
  # Genome labels (italic, separated)
  geom_text_repel(
    data = scores_labels,
    aes(PC1, PC2, label = Species),
    size = 3.5,
    fontface = "italic",
    max.overlaps = 50,
    box.padding = 0.7,
    point.padding = 0.8,
    force = 2,
    force_pull = 0.3,
    segment.color = "grey65",
    segment.size = 0.25, show.legend = F
  ) +
  
  # Category vectors
  geom_segment(
    data = loadings_plot,
    aes(
      x = 0, y = 0,
      xend = PC1, yend = PC2,
      color = Lifestyle
    ),
    arrow = arrow(length = unit(0.25, "cm")),
    linewidth = 0.7
  ) +
  
  # Category labels (NO padding)
  geom_text_repel(
    data = loadings_plot,
    aes(PC1, PC2, label = Category, color = Lifestyle),
    size = 3.5,
    max.overlaps = 100,
    box.padding = 0, show.legend = F,
    point.padding = 0
  ) +
  
  scale_color_manual(
    values = lifestyle_colors,
    breaks = c("Pathogen", "Symbiont", "Both"),
    labels = c("Pathogen", "Symbiont", "Both"),
    name = "Strategy"
  ) +
  
  theme_bw(base_size = 12) +
  labs(
    title = "PCA of correlated genomes based on KEGG repertoire",
    x = paste0(
      "PC1 (",
      round(100 * summary(pca)$importance[2,1], 1),
      "%)"
    ),
    y = paste0(
      "PC2 (",
      round(100 * summary(pca)$importance[2,2], 1),
      "%)"
    )
  ) +
  
  theme(
    panel.grid = element_blank(),
    legend.position = "right",
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14)
  )
pca.corr.kegg

tiff("KEGG_PCA_correlated_genomes_minimal_2026.tiff", units="cm", width=16, height=12, res = 300, compression = 'lzw')
pca.corr.kegg
dev.off()



## Indicator taxon plots
# Plot frequencies of pathways - KEGG + Indicator
df_heatmap_indic <- indic.gene.dbb.kegg.target %>%
  group_by(Species, Category, Correlation, Lifestyle) %>%
  summarise(n = n(), .groups = "drop")

# Order Species alphabetically within each Correlation (Indicator heatmaps)
species_order_indic <- indic.gene.dbb.kegg.target %>%
  select(Species, Correlation) %>%
  distinct() %>%
  arrange(Correlation, Species) %>%
  pull(Species)

# Order Categories alphabetically (reversed) within each Lifestyle (y axis)
category_order_indic <- indic.gene.dbb.kegg.target %>%
  select(Category, Lifestyle) %>%
  distinct() %>%
  group_by(Lifestyle) %>%
  arrange(Lifestyle, desc(Category)) %>%  # reverse alphabetical *within* Lifestyle
  ungroup() %>%
  pull(Category)

df_heatmap_indic <- df_heatmap_indic %>%
  mutate(
    Species  = factor(Species, levels = unique(species_order_indic)),
    Category = factor(Category, levels = unique(category_order_indic))
  )

category_label_df_indic <- indic.gene.dbb.kegg.target %>%
  select(Category, Lifestyle) %>%
  distinct() %>%
  mutate(
    label_colored = paste0(
      "<span style='color:",
      lifestyle_colors[Lifestyle],
      "'>",
      Category,
      "</span>"
    )
  )

f_heatmap_indic <- ggplot(df_heatmap_indic, aes(x = Species, y = Category, fill = n)) +
  geom_tile(colour = "white", linewidth = 0.3) +
  scale_fill_gradient(
    low = "#EAF3DE", high = "#173404",
    name = "KEGG count", na.value = "grey95"
  ) +
  scale_y_discrete(
    drop = FALSE,
    labels = setNames(
      category_label_df_indic$label_colored,
      category_label_df_indic$Category
    )
  ) +
  facet_grid(. ~ Correlation, scales = "free_x", space = "free_x") +
  labs(
    x = NULL, y = NULL,
    title = "KEGG repertoire per indicator genome"
  ) +
  theme_bw(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
    axis.text.y = element_markdown(size = 12),
    strip.text = element_text(face = "bold", size = 7),
    strip.background = element_rect(fill = "grey92"),
    legend.position = "right",
    legend.key.size = unit(0.5, "cm"),
    panel.grid = element_blank(),
    plot.title = element_text(size = 14, face = "bold")
  )

f_heatmap_indic

pdf("KEGG_heatmap_indicator_genomes_2026.pdf", width = 12, height = 6.5)
f_heatmap_indic
dev.off()


# Order Species by number of Categories present
species_order_by_categories_indic <- df_heatmap_indic %>%
  distinct(Species, Category) %>%
  count(Species, name = "n_categories") %>%
  arrange(desc(n_categories)) %>%
  pull(Species)

df_heatmap_indic <- df_heatmap_indic %>%
  mutate(Species = factor(Species, levels = species_order_by_categories_indic))

n_heatmap_indic <- ggplot(df_heatmap_indic, aes(x = Species, y = Category, fill = n)) +
  geom_tile(colour = "white", linewidth = 0.3) +
  scale_fill_gradient(
    low = "#EAF3DE", high = "#173404",
    name = "KEGG count", na.value = "grey95"
  ) +
  scale_y_discrete(
    drop = FALSE,
    labels = setNames(
      category_label_df_indic$label_colored,
      category_label_df_indic$Category
    )
  ) +
  labs(
    x = NULL, y = NULL,
    title = "KEGG repertoire per indicator genome"
  ) +
  theme_bw(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14, face = "italic"),
    axis.text.y = element_markdown(size = 12),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 10),
    legend.position = "right",
    legend.key.size = unit(0.5, "cm"),
    panel.grid = element_blank(),
    plot.title = element_text(size = 14, face = "bold")
  )

n_heatmap_indic

pdf("KEGG_heatmap_indicator_genomes_2026.pdf", width = 12, height = 6.5)

n_heatmap_indic
dev.off()


## Plot as ordination
df_heatmap_indic <- as.data.frame(df_heatmap_indic)

kegg_mat_indic <- df_heatmap_indic %>%
  group_by(Species, Category) %>%
  summarise(n = sum(n), .groups = "drop") %>%
  tidyr::pivot_wider(
    names_from  = Category,
    values_from = n,
    values_fill = list(n = 0)
  ) %>%
  tibble::column_to_rownames("Species")

pca_indic <- prcomp(kegg_mat_indic, center = TRUE, scale. = TRUE)
summary(pca_indic)$importance[2, 1:2]

scores_indic <- as.data.frame(pca_indic$x[, 1:2])
scores_indic$Species <- rownames(scores_indic)

loadings_indic <- as.data.frame(pca_indic$rotation[, 1:2])
loadings_indic$Category <- rownames(loadings_indic)

arrow_scale <- 5
loadings_plot_indic <- loadings_indic %>%
  mutate(
    PC1 = PC1 * arrow_scale,
    PC2 = PC2 * arrow_scale
  ) %>%
  filter(abs(PC1) + abs(PC2) > 2)

category_lifestyle_indic <- df_heatmap_indic %>%
  select(Category, Lifestyle) %>%
  distinct()
loadings_plot_indic <- loadings_plot_indic %>%
  left_join(category_lifestyle_indic, by = "Category")

scores_labels_indic <- scores_indic %>%
  filter(PC1 < -2.5)

pca.indic.kegg <- ggplot() +
  
  # Genome points
  geom_point(
    data = scores_indic,
    aes(PC1, PC2),
    size = 3.5,
    color = "black"
  ) +
  
  # Genome labels (italic, separated)
  geom_text_repel(
    data = scores_labels_indic,
    aes(PC1, PC2, label = Species),
    size = 3.5,
    fontface = "italic",
    max.overlaps = 50,
    box.padding = 0.8,     # ??? more space around labels
    point.padding = 0.8,   # ??? more space from points
    force = 3,             # ??? stronger repulsion between labels
    force_pull = 0.25,     # ??? slightly weaker attraction to points
    segment.color = "grey65",
    segment.size = 0.25,
    show.legend = FALSE
  ) +
  
  # Category vectors
  geom_segment(
    data = loadings_plot_indic,
    aes(
      x = 0, y = 0,
      xend = PC1, yend = PC2,
      color = Lifestyle
    ),
    arrow = arrow(length = unit(0.25, "cm")),
    linewidth = 0.7
  ) +
  
  # Category labels (no padding)
  geom_text_repel(
    data = loadings_plot_indic,
    aes(PC1, PC2, label = Category, color = Lifestyle),
    size = 3.5,
    max.overlaps = 100,
    box.padding = 0,
    point.padding = 0,
    show.legend = FALSE
  ) +
  
  scale_color_manual(
    values = lifestyle_colors,
    breaks = c("Pathogen", "Symbiont", "Both"),
    labels = c("Pathogen", "Symbiont", "Both"),
    name = "Strategy"
  ) +
  
  theme_bw(base_size = 12) +
  labs(
    title = "PCA of indicator genomes based on KEGG repertoire",
    x = paste0(
      "PC1 (",
      round(100 * summary(pca_indic)$importance[2,1], 1),
      "%)"
    ),
    y = paste0(
      "PC2 (",
      round(100 * summary(pca_indic)$importance[2,2], 1),
      "%)"
    )
  ) +
  
  theme(
    panel.grid = element_blank(),
    legend.position = "right",
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14)
  )

pca.indic.kegg

tiff("KEGG_PCA_indicator_genomes_minimal_2026.tiff", units = "cm",
  width = 16, height = 12, res = 300, compression = "lzw")
pca.indic.kegg
dev.off()
