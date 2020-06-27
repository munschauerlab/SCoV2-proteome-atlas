library(STRINGdb)
library(data.table)
library(dplyr)

# Import database
string_db <- STRINGdb$new( version="10", species=9606, score_threshold=50, input_directory="../string")

# Import gene names
rap <- fread("../data/MM_proteomics_data_limmaFDR.txt") %>%
  filter(adj.P.Val.SARS.CoV.over.RMRP < 0.2)

# Map to STRING Ids
rap_mapped <- string_db$map( rap, "geneSymbol", removeUnmappedRows = TRUE ) 

# See what fell out
rap$geneSymbol[!(rap$geneSymbol %in% rap_mapped$geneSymbol)]

# Build interaction network
query_proteins <- rap_mapped$STRING_id
trans_vec <- c(rap_mapped$geneSymbol); names(trans_vec) <- query_proteins
edges <- string_db$get_interactions(unique(query_proteins))[,c("from", "to", "combined_score")] %>%
  arrange(desc(combined_score))
colnames(edges) <- c("A", "B", "score")
edges$node1 <- trans_vec[as.character(edges$A)]
edges$node2 <- trans_vec[as.character(edges$B)]

# QC genes that may be missing
qc_genes <- c("SCFD1", "RAB6C", "RAB7A")

edges %>% filter(node1 %in% qc_genes | node2 %in% qc_genes)

# Visualize network
levels = unique(c(edges[["node1"]], edges[["node2"]]))
edges$node1 <- factor(edges$node1, levels = levels)
edges$node2 <- factor(edges$node2, levels = levels)

sm <- sparseMatrix(i = c(as.numeric(edges[["node1"]]), length(levels)), 
                   j = c(as.numeric(edges[["node2"]]), length(levels)), 
                   x  = c(edges[["score"]], 0)) %>% data.matrix()
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

set.seed(3)
pAI <- ggnet2(net, label =TRUE, color = "color",  legend.position = "bottom", 
              palette = c("Other" = "grey", "CNBP" = "red", "RP" =  "dodgerblue"), size = 3, label.size = 4)
cowplot::ggsave2(pAI, file = "../output/string/String_full_network_27June2020.pdf", width = 5, height = 5)

