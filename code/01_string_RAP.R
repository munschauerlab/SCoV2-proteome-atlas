library(data.table)
library(dplyr)

# Import database
if(!exists("string_v11")){
  source('00_string_v11.R')
}

# Import gene names
rap <- fread("../data/MM_proteomics_data_limmaFDR.txt") %>%
  filter(adj.P.Val.SARS.CoV.over.RMRP < 0.2 & geneSymbol %ni% viral_genes)
rap_genes <- rap %>% pull(geneSymbol)

# Make string network
edges <- pull_string_v11(min_score = 750, genes_keep = rap_genes)

# QC genes that may be missing
qc_genes <- c("SCFD1", "RAB6C", "RAB7A")

edges %>% filter(node1 %in% qc_genes | node2 %in% qc_genes)

# Visualize network
levels = unique(c(edges[["node1"]], edges[["node2"]]))
edges$node1 <- factor(edges$node1, levels = levels)
edges$node2 <- factor(edges$node2, levels = levels)

sm <- sparseMatrix(i = c(as.numeric(edges[["node1"]]), length(levels)), 
                   j = c(as.numeric(edges[["node2"]]), length(levels)), 
                   x  = c(edges[["combined_score"]], 0)) %>% data.matrix()
rownames(sm) <- levels
colnames(sm) <- levels

net <- network(sm)

# Show the full network
color_vec <- case_when(
  grepl("^RP", levels) ~ "RP", 
  levels == "CNBP" ~ "CNBP",
  TRUE ~ "Other"
)
net %v% "color" = color_vec

i = 1
set.seed(i)
pAI <- ggnet2(net, label = TRUE, color = "color",  legend.position = "bottom", mode = "fruchtermanreingold",
              layout.par = list(repulse.rad = 2000),
              palette = c("Other" = "grey", "CNBP" = "red", "RP" =  "dodgerblue"), size = 3, label.size = 3)
cowplot::ggsave2(pAI, file = paste0("../output/string/string_seed/String_RAPnetwork",as.character(i),".pdf"), width = 5, height = 5)
  

