library(GGally)
library(network)
library(data.table)
library(Matrix)
library(dplyr)

mm <- fread("../data/MM_proteomics_data_limmaFDR.txt")
mm_fdr20 <- mm %>% filter(adj.P.Val.SARS.CoV.over.RMRP < 0.2) %>% pull(geneSymbol)
edges <- fread("../data/genets_fdr20_network/DownloadEdges_424468.csv")
edges <- edges %>% filter(score > 0.3)
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
rb <- c("RPS2", "RPS14", "RPL15", "RPS10","NUDT3", "RPL6", "RPS3", "RPL13", "RPS4X", "RPS26", "RPS11", "RPL21", "RPS12", "RPS5", "RPL8", "RPL18A", "RPL3", "RPL36A", "RPL28", "RPL7A")
rtv <- c("C19orf66", "DDX1", "LSM14A", "DDX3X", "PCBP2", "ACTA2", "CFL1")
vp <- c("RAB6A", "RAB7B", "RAB1A", "EGFR", "PPIA", "SND1", "APOE", "HNRNPA1", "PCBP2", "DDX3X", "SYNCRIP", "EIF4G1", "EIF3G", "EIF4H", "EIF3H")
iep <- c("RAB6A", "RAB7B", "STOM", "C19orf66", "LSM14A", "PCBP2", "DDX3X", "ACTR2", "ACTB", "PPIA", "EEF2", "EEF1A1", "GPI")

# Annotate colors
net %v% "rb" = (levels %in% rb)
net %v% "rtv" = (levels %in% rtv)
net %v% "vp" = (levels %in% vp)
net %v% "iep" = (levels %in% iep)

set.seed(3)
p_rb <- ggnet2(net, label = rb, color = "rb",  legend.position = "none", 
       palette = c("FALSE" = "grey", "TRUE" = "red"), size = 3, label.size = 4)
cowplot::ggsave2(p_rb, file = "../output/genets/network_rb.pdf", width = 5, height = 5)

set.seed(3)
p_rtv <- ggnet2(net, label = rtv, color = "rtv",  legend.position = "none", 
               palette = c("FALSE" = "grey", "TRUE" = "red"), size = 3, label.size = 4)
cowplot::ggsave2(p_rtv, file = "../output/genets/network_rtv.pdf", width = 5, height = 5)

set.seed(3)
p_vp <- ggnet2(net, label = vp, color = "vp",  legend.position = "none", 
               palette = c("FALSE" = "grey", "TRUE" = "red"), size = 3, label.size = 4)
cowplot::ggsave2(p_vp, file = "../output/genets/network_vb.pdf", width = 5, height = 5)

set.seed(3)
p_iep <- ggnet2(net, label = iep, color = "iep",  legend.position = "none", 
               palette = c("FALSE" = "grey", "TRUE" = "red"), size = 3, label.size = 4)
cowplot::ggsave2(p_iep, file = "../output/genets/network_iep.pdf", width = 5, height = 5)


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
cowplot::ggsave2(pAI, file = "../output/genets/full_network_27June2020.pdf", width = 5, height = 5)

