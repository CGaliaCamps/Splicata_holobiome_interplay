rm(list = ls())

# Indicator ASV visualisation
# Stacked bar: taxonomic composition per indicator group, host-association highlighted

library(tidyverse)
library(patchwork)
library(scales)

setwd("E:/Splicata/16S-RNA-DNA/Micro_function")

df_raw <- read.csv("Host_association_indicator_ASVs.csv",
                   stringsAsFactors = FALSE, na.strings = c("", "NA"))


df_raw %>%
  count(Host_association) %>%
  mutate(proportion = n / sum(n)) # 2.9% certain (Francisellaceae and Nitrosopumilus), 48.6% certain, 32.9% unlikely, 12.7 no association

#0.02857143 certain + 0.48571429 likely = 51.4% certain + likely (host-associated)

assoc_levels  <- c("Certain", "Likely", "Unlikely", "No association")
assoc_colours <- c(
  "Certain"        = "#1D6A47",
  "Likely"         = "#5AAE8C",
  "Unlikely"       = "#D4956A",
  "No association" = "#BBBBBB"
)

class_colours <- c(
  "Alphaproteobacteria" = "#2E5FA8",
  "Gammaproteobacteria" = "#8BB0E8",
  "Bacteroidia"         = "#ED7D31",
  "Cyanobacteriia"      = "#70AD47",
  "Actinobacteria"      = "#9E5FAB",
  "Acidimicrobiia"      = "#C89FCC",
  "Coriobacteriia"      = "#E0C0E8",
  "Nitrososphaeria"     = "#E8C53A",   # Crenarchaeota
  "Polyangia"           = "#37A9CC",   # Myxococcota
  "OM190"               = "#2EAAB8",   # Planctomycetota
  "Nitrospiria"         = "#3A8C6E",   # Nitrospirota
  "Gracilibacteria"     = "#E84040",   # Patescibacteria
  "Verrucomicrobiae"    = "#8A7A5A",
  "Unclassified"        = "#CCCCCC"
)

# Keeping codes but adding readable labels
indicator_labels <- c(
  "BLA75"  = "BLA75",
  "VIL75"  = "VIL75",
  "BAR75"  = "BAR75",
  "2AA75"  = "2AA75",
  "2AB75"  = "2AB75",
  "4AB75"  = "4AB75",
  "11AB75" = "11AB75",
  "11BB75" = "11BB75",
  "16AA75" = "16AA75",
  "16AB75" = "16AB75"
)

# Order: put BLA75 first (most ASVs), then others by n desc
indicator_order <- names(sort(table(df_raw$indicator), decreasing = TRUE))

df <- df_raw %>%
  mutate(
    Class            = if_else(is.na(Class), "Unclassified", Class),
    Family           = if_else(is.na(Family), "Unclassified", Family),
    Host_association = factor(Host_association, levels = assoc_levels),
    indicator        = factor(indicator, levels = indicator_order),
    Class            = factor(Class, levels = names(class_colours))
  )

bar_df <- df %>%
  count(indicator, Class, Host_association) %>%
  mutate(Class = factor(Class, levels = names(class_colours)))

assoc_strip <- df %>%
  count(indicator, Host_association) %>%
  group_by(indicator) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

p_top <- ggplot(assoc_strip,
                aes(x = indicator, y = prop, fill = Host_association)) +
  geom_col(width = 0.72, colour = "white", linewidth = 0.3) +
  scale_fill_manual(values = assoc_colours, name = "Host\nassociation",
                    guide = guide_legend(ncol = 1,
                                         keyheight = unit(0.85, "lines"))) +
  scale_x_discrete(labels = indicator_labels) +
  scale_y_continuous(labels = percent_format(accuracy = 1),
                     expand = expansion(mult = c(0, 0.12)),
                     breaks = c(0, 0.5, 1)) +
  labs(x = NULL, y = "Proportion") +
  theme_classic(base_size = 11) +
  theme(
    axis.text.x       = element_blank(),
    axis.text.y       = element_text(size = 10),
    axis.ticks.x      = element_blank(),
    legend.position   = "right",
    legend.text       = element_text(size = 10),
    legend.title      = element_text(size = 10, face = "bold"),
    panel.grid.major.y = element_line(colour = "grey90", linewidth = 0.3),
    plot.margin       = margin(6, 6, 0, 6)
  )

p_bot <- ggplot(bar_df, aes(x = indicator, y = n, fill = Class)) +
  geom_col(width = 0.72, colour = "white", linewidth = 0.3) +
  scale_fill_manual(values = class_colours, name = "Class",
                    guide = guide_legend(ncol = 1,
                                         keyheight = unit(0.82, "lines"))) +
  scale_x_discrete(labels = indicator_labels) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.04)),
                     breaks = pretty_breaks(n = 6)) +
  labs(
    x = "Indicator group",
    y = "Number of ASVs"
  ) +
  theme_classic(base_size = 11) +
  theme(
    axis.text.y       = element_text(size = 10),
    axis.text.x       = element_text(size = 10, angle = 30, hjust = 1),
    legend.position   = "right",
    legend.text       = element_text(size = 10),
    legend.title      = element_text(size = 10, face = "bold"),
    panel.grid.major.y = element_line(colour = "grey90", linewidth = 0.3),
    plot.margin       = margin(0, 6, 6, 6)
  )

figx_indicator <- p_top / p_bot +
  plot_layout(heights = c(1, 2.8), guides = "collect") +
  plot_annotation(
    theme    = theme(
      plot.title    = element_text(face = "bold", size = 12),
      plot.subtitle = element_text(size = 9, colour = "grey40")
    )
  ) &
  theme(legend.position = "right")

figx_indicator

ggsave("Taxonomy_by_indicator_2026_legend.pdf",
       plot        = figx_indicator,
       width       = 7, height = 5)

ggsave("Taxonomy_by_indicator_2026.pdf",
       plot        = figx_indicator,
       width       = 6.5, height = 3.75)


