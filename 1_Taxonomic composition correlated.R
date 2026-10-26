rm(list = ls())

library(tidyverse)
library(ggalluvial)    
library(patchwork)     
library(scales)

# 1. Data
setwd("E:/Splicata/16S-RNA-DNA/Micro_function")

df_raw <- read.csv("Host_association_correlated_ASVs.csv",
                   stringsAsFactors = FALSE, na.strings = c("", "NA"))
head(df_raw)

df_raw %>%
  count(Host_association) %>%
  mutate(proportion = n / sum(n)) # 3.6% certain (Francisellaceae and Nitrosopumilus), 54.5% certain, 27.3% unlikely, 14.5% no association

#0.03636364 certain + 0.54545455 likely = 58.2% certain + likely (host-associated) 

# 2. Factor levels and palettes
assoc_levels <- c("Certain", "Likely", "Unlikely", "No association")
assoc_colours <- c(
  "Certain"        = "#1D6A47",   # deep teal-green
  "Likely"         = "#5AAE8C",   # mid teal
  "Unlikely"       = "#D4956A",   # warm amber-orange
  "No association" = "#BBBBBB"    # neutral grey
)

# Correlation groups - rename for readable labels
corr_labels <- c(
  "pop"           = "Population",
  "chr11"         = "Chromosome 11",
  "chr11_onlychr" = "Chromosome 11 only",
  "chr2"          = "Chromosome 2",
  "chr2_onlychr"  = "Chromosome 2 only",
  "chr4"          = "Chromosome 4"
)

corr_order <- c("pop", "chr11", "chr11_onlychr", "chr2", "chr2_onlychr", "chr4")

df <- df_raw %>%
  mutate(
    Host_association = factor(Host_association, levels = assoc_levels),
    Correlation      = factor(Correlation, levels = corr_order),
    # Combine Phylum:Class for the stacked bar fill
    Phylum_Class     = paste0(Phylum, "\n(", Class, ")"),
    # Tidy family label
    Family           = if_else(is.na(Family), "Unclassified", Family)
  )

df.chronly = df[df$Correlation %in% c("pop", "chr11_onlychr", "chr2_onlychr"),]
df.all = df[df$Correlation %in% c("pop", "chr11", "chr2", "chr4"),]

# 3. Phylum colour palette
phylum_colours <- c(
  "Proteobacteria"    = "#4472C4",
  "Bacteroidota"      = "#ED7D31",
  "Cyanobacteria"     = "#70AD47",
  "Actinobacteriota"  = "#9E5FAB",
  "Crenarchaeota"     = "#E8C53A",
  "Myxococcota"       = "#37A9CC",
  "Patescibacteria"   = "#E84040",
  "Verrucomicrobiota" = "#8A7A5A",
  "Unassigned"        = "#CCCCCC"
)

# Class within Proteobacteria - distinguish Alpha vs. Gamma with hue shift
class_colours <- c(
  # Proteobacteria
  "Alphaproteobacteria" = "#2E5FA8",
  "Gammaproteobacteria" = "#8BB0E8",
  # Bacteroidota
  "Bacteroidia"         = "#ED7D31",
  # Cyanobacteria
  "Cyanobacteriia"      = "#70AD47",
  # Actinobacteriota
  "Actinobacteria"      = "#9E5FAB",
  "Acidimicrobiia"      = "#C89FCC",
  # Others
  "Nitrososphaeria"     = "#E8C53A",
  "Polyangia"           = "#37A9CC",
  "Gracilibacteria"     = "#E84040",
  "Verrucomicrobiae"    = "#8A7A5A",
  "Unassigned"          = "#CCCCCC"
)

# 4. Stacked bar: taxonomy ? correlation ? association
# Strategy: fill = Class (taxonomy), alpha encodes association confidence
# Association shown via a second thin bar / annotation strip on top

# Summarise counts
bar_df <- df.chronly %>%
  count(Correlation, Class, Host_association) %>%
  mutate(
    Correlation = factor(Correlation, levels = corr_order),
    Class       = factor(Class, levels = names(class_colours))
  )

# Association summary per correlation (for annotation strip)
assoc_strip <- df.chronly %>%
  count(Correlation, Host_association) %>%
  group_by(Correlation) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup() %>%
  mutate(Correlation = factor(Correlation, levels = corr_order))

# Main stacked bar
p1_main <- ggplot(bar_df, aes(x = Correlation, y = n, fill = Class)) +
  geom_col(width = 0.7, colour = "white", linewidth = 0.3) +
  scale_fill_manual(values = class_colours, name = "Class",
                    guide = guide_legend(ncol = 1, keyheight = unit(0.85, "lines"))) +
  scale_x_discrete(labels = corr_labels) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(
    x = "Correlation group",
    y = "Number of ASVs",
  ) +
  theme_classic(base_size = 11) +
  theme(
    legend.position   = "right",
    legend.text       = element_text(size = 10),
    legend.title      = element_text(size = 10, face = "bold"),
    axis.text.x       = element_text(size = 10, angle = 45, vjust = 1, hjust = 1, lineheight = 1.1),
    plot.title        = element_text(face = "bold", size = 12),
    plot.subtitle     = element_text(size = 9, colour = "grey40"),
    panel.grid.major.y = element_line(colour = "grey90", linewidth = 0.3)
  )

# Association strip below (proportional stacked bar, same x-axis)
p1_strip <- ggplot(assoc_strip, aes(x = Correlation, y = prop, fill = Host_association)) +
  geom_col(width = 0.7, colour = "white", linewidth = 0.3) +
  scale_fill_manual(values = assoc_colours, name = "Host association",
                    guide = guide_legend(ncol = 1, keyheight = unit(0.85, "lines"))) +
  scale_x_discrete(labels = corr_labels) +
  scale_y_continuous(labels = percent_format(), expand = expansion(mult = c(0, 0.02))) +
  labs(x = NULL, y = "Proportion") +
  theme_classic(base_size = 11) +
  theme(
    legend.position   = "right",
    legend.text       = element_text(size = 10),
    legend.title      = element_text(size = 10, face = "bold"),
    axis.text.x       = element_blank(),
    axis.ticks.x      = element_blank(),
    panel.grid.major.y = element_line(colour = "grey90", linewidth = 0.3)
  )

figx <- p1_strip / p1_main +
  plot_layout(heights = c(1, 3), guides = "collect") &
  theme(legend.position = "right")

figx

ggsave("Taxonomy_by_correlated_onlychr_2026.pdf",
       plot        = figx,
       width       = 6.5, height = 3.75)




## All included
# Summarise counts
bar_df <- df.all %>%
  count(Correlation, Class, Host_association) %>%
  mutate(
    Correlation = factor(Correlation, levels = corr_order),
    Class       = factor(Class, levels = names(class_colours))
  )

# Association summary per correlation (for annotation strip)
assoc_strip <- df.all %>%
  count(Correlation, Host_association) %>%
  group_by(Correlation) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup() %>%
  mutate(Correlation = factor(Correlation, levels = corr_order))

# Main stacked bar
p1_main <- ggplot(bar_df, aes(x = Correlation, y = n, fill = Class)) +
  geom_col(width = 0.7, colour = "white", linewidth = 0.3) +
  scale_fill_manual(values = class_colours, name = "Class",
                    guide = guide_legend(ncol = 1, keyheight = unit(0.85, "lines"))) +
  scale_x_discrete(labels = corr_labels) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(
    x = "Correlation group",
    y = "Number of ASVs",
  ) +
  theme_classic(base_size = 11) +
  theme(
    legend.position   = "right",
    legend.text       = element_text(size = 10),
    legend.title      = element_text(size = 10, face = "bold"),
    axis.text.x       = element_text(size = 10, angle = 45, vjust = 1, hjust = 1, lineheight = 1.1),
    plot.title        = element_text(face = "bold", size = 12),
    plot.subtitle     = element_text(size = 9, colour = "grey40"),
    panel.grid.major.y = element_line(colour = "grey90", linewidth = 0.3)
  )

# Association strip below (proportional stacked bar, same x-axis)
p1_strip <- ggplot(assoc_strip, aes(x = Correlation, y = prop, fill = Host_association)) +
  geom_col(width = 0.7, colour = "white", linewidth = 0.3) +
  scale_fill_manual(values = assoc_colours, name = "Host association",
                    guide = guide_legend(ncol = 1, keyheight = unit(0.85, "lines"))) +
  scale_x_discrete(labels = corr_labels) +
  scale_y_continuous(labels = percent_format(), expand = expansion(mult = c(0, 0.02))) +
  labs(x = NULL, y = "Proportion") +
  theme_classic(base_size = 11) +
  theme(
    legend.position   = "right",
    legend.text       = element_text(size = 10),
    legend.title      = element_text(size = 10, face = "bold"),
    axis.text.x       = element_blank(),
    axis.ticks.x      = element_blank(),
    panel.grid.major.y = element_line(colour = "grey90", linewidth = 0.3)
  )

figx <- p1_strip / p1_main +
  plot_layout(heights = c(1, 3), guides = "collect") &
  theme(legend.position = "right")

figx

ggsave("Taxonomy_by_correlated_2026_legend.pdf",
       plot        = figx,
       width       = 7, height = 5)

ggsave("Taxonomy_by_correlated_2026.pdf",
       plot        = figx,
       width       = 6.5, height = 3.75)
