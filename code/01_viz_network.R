library(GGally)
library(network)
library(data.table)
library(Matrix)
library(dplyr)

mm <- fread("../data/MM_proteomics_data_limmaFDR.txt")
mm_fdr20 <- mm %>% filter(adj.P.Val.SARS.CoV.over.RMRP < 0.2) %>% pull(geneSymbol)
edges <- fread("../data/genets_fdr20_network/DownloadEdges_424468.csv")
edges <- edges %>% filter(score > 0.4)
edges <- edges[edges[["node1"]] %in% mm_fdr20 & edges[["node2"]] %in% mm_fdr20,]


levels = unique(c(edges[["node1"]], edges[["node2"]]))
edges$node1 <- factor(edges$node1, levels = levels)
edges$node2 <- factor(edges$node2, levels = levels)

sm <- sparseMatrix(i = c(as.numeric(edges[["node1"]]), length(levels)), 
                   j = c(as.numeric(edges[["node2"]]), length(levels)), 
                   x  = c(edges[["score"]], 0)) %>% data.matrix()
rownames(sm) <- levels
colnames(sm) <- levels

net <- network(sm)
color_vec <- case_when(
  grepl("^RP", levels) ~ "RP", 
  levels == "CNBP" ~ "CNBP",
  TRUE ~ "Other"
)
net %v% "color" = color_vec

set.seed(3)
pAI <- ggnet2(net, label =TRUE, color = "color",  legend.position = "bottom", 
       palette = c("Other" = "grey", "CNBP" = "red", "RP" =  "dodgerblue"), size = 3, label.size = 4)
cowplot::ggsave2(pAI, file = "../output/toy_network.pdf", width = 5, height = 5)
