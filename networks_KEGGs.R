
library(stringr)
library(dplyr)
library(tidyr)
library(purrr)

setwd("E://Splicata/16S-RNA-DNA/networks")


txt <- readLines("host-symbiont_KEGG_module.txt", warn = FALSE)
txt <- paste(txt, collapse = "\n")

# split by module
modules <- str_split(txt, "\n(?=M\\d{5})")[[1]]

parse_module <- function(block){
  
  # module ID
  module <- str_extract(block, "M\\d{5}")
  
  # pathway name (line after module line)
  pathway <- str_match(block, paste0(module, "\\s+([^\\(\\n]+)"))[,2]
  pathway <- str_trim(pathway)
  
  # split by K entries
  k_blocks <- str_split(block, "\nK\\d{5}")[[1]]
  
  # remove first chunk (header before first K)
  k_blocks <- k_blocks[-1]
  
  k_ids <- str_extract_all(block, "K\\d{5}")[[1]]
  
  map2_dfr(k_blocks, k_ids, function(kb, kid){
    
    # ASVs
    asvs <- str_match_all(kb, "corr_iASV:([^\\n]+)")[[1]][,2]
    asvs <- if(length(asvs) > 0)
      str_split(asvs, ",\\s*")[[1]] else character(0)
    
    asvs <- str_trim(asvs)
    
    # DEGs (SPLI)
    degs <- str_match_all(kb, "host_Splicata:([^\\n]+)")[[1]][,2]
    degs <- degs[str_detect(degs, "SPLI_")]
    degs <- if(length(degs) > 0)
      str_split(degs, ",\\s*")[[1]] else character(0)
    
    degs <- str_trim(degs)
    
    tibble(
      module = module,
      pathway = pathway,
      KEGG = kid
    ) %>%
      bind_rows(
        tibble(module = module, pathway = pathway, KEGG = kid,
               type = "ASV", ID = asvs),
        tibble(module = module, pathway = pathway, KEGG = kid,
               type = "DEG", ID = degs)
      ) %>%
      filter(!is.na(ID))
  })
}

df <- map_dfr(modules, parse_module)

write.table(df, "KEGG-ASV_Modules.txt", sep ="\t")


edges_kegg <- df %>%
  filter(type %in% c("ASV", "DEG")) %>%
  group_by(module, KEGG) %>%
  summarise(
    ASVs = list(unique(ID[type == "ASV"])),
    DEGs = list(unique(ID[type == "DEG"])),
    .groups = "drop"
  ) %>%
  rowwise() %>%
  mutate(edges = list(
    if (length(ASVs) > 0 && length(DEGs) > 0)
      expand.grid(from = ASVs, to = DEGs)
    else NULL
  )) %>%
  tidyr::unnest(edges) %>%
  ungroup() %>%
  select(from, to, KEGG, module)

g <- graph_from_data_frame(edges_kegg, directed = FALSE)
kegg_levels <- unique(E(g)$KEGG)
module_levels <- unique(E(g)$module)
V(g)$type <- ifelse(grepl("^ASV", V(g)$name), "ASV", "DEG")

kegg_cols <- setNames(
  c("#FF7F00", "#7FFF00","#007FFF","#00FF7F"),
  kegg_levels
)

E(g)$color <- kegg_cols[E(g)$KEGG]
E(g)$color <- adjustcolor(E(g)$color, alpha.f = 0.5)

svg("network_KEGG.svg")
plot(
  g,
  vertex.size = 10,
  vertex.color = ifelse(V(g)$type == "ASV", "grey", "white"),
  vertex.label.cex = NA,
  #vertex.label.color="black",
  edge.color = E(g)$color,
  edge.width = 2
)

legend(
  "topleft",
  legend = names(kegg_cols),
  col = kegg_cols,
  lty = 1,
  cex = 1
)
dev.off()


module_cols <- setNames(
  c("#FF0000","yellow2", "#00FF00","#7F00FF", "#00FFFF"),
  module_levels
)

E(g)$color <- module_cols[E(g)$module]
E(g)$color <- adjustcolor(E(g)$color, alpha.f = 0.5)

svg("network_module.svg")
plot(
  g,
  vertex.size = 10,
  vertex.color = ifelse(V(g)$type == "ASV", "grey", "white"),
  vertex.label.cex = 1,
  vertex.label.color="black",
  edge.color = E(g)$color,
  edge.width = 2
)

legend(
  "topleft",
  legend = names(module_cols),
  col = module_cols,
  lty = 1,
  cex = 1
)
dev.off()

