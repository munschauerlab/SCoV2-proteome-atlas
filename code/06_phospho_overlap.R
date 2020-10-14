library(data.table)
library(dplyr)
library(Matrix)
library(GGally)
library(network)

# Import database
if(!exists("string_v11")){
  source('00_string_v11.R')
}

# Import gene names
rap <- fread("../data/SCoV2_RAPms.txt") %>%
  filter(adj.P.Val.SCoV2.over.RMRP < 0.2 & species == "HOMO SAPIENS" & logFC.SCoV2.over.RMRP > 0)
rap_genes <- rap %>% pull(geneSymbol)

# Make string network
edges <- pull_string_v11(min_score = 550, genes_keep = rap_genes)
write.xlsx(edges, "../output/TableS4_interactome_only_network.xlsx",
           sheetName = "Interactome", 
           col.names = TRUE, row.names = FALSE, append = TRUE)

# QC genes that may be missing
qc_genes <- c("SCFD1", "RAB6C", "RAB7A", "GDI2")

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
vesicle <- c("RAB7A","RAB6A","RAB1A","SCFD1","USO1", "GDI2")
viral_rna <- c("ATP1A1","CAPRIN1","CFL1","CSDE1","DDX1","DDX3X","EEF1A1","EEF2","EIF3E","EIF3G","EIF3L","EIF4B","EIF4G1","EIF4H","EIF5A","G3BP1","G3BP2",
               "HNRNPA1","HNRNPA2B1","IGF2BP1","IGF2BP2","LIN28B","LSM14A","MOV10","PABPC1","PCBP2","PEBP1","PFN1","PPIA","RPL13","RPL15","RPL18A","RPL21",
               "RPL28","RPL3","RPL36A","RPL6","RPL7A","RPL8","RPS11","RPS12","RPS14","RPS2","RPS26","RPS3","RPS4X","RPS5","SND1","SYNCRIP","UPF1","YBX1")

color_vec <- case_when(
  grepl("^RP", levels) ~ "RP",
  levels %in% viral_rna ~ "viral_rna_bind",
  TRUE ~ "Other"
)
table(color_vec)

net %v% "color" = color_vec

# Now do size
n_size_vec = rap$logFC.SCoV2.over.RMRP; names(n_size_vec) <- as.character(rap$geneSymbol)
size_vec <- unname(n_size_vec[levels])
net %v% "size_of_node" = size_vec

Bouhaddou <- readxl::read_xlsx("../data/Bouhaddou2020_TS1.xlsx") %>%
  data.frame() %>%
  filter(Inf_00Hr.adj.pvalue < 0.05 | Inf_02Hr.adj.pvalue < 0.05 | Inf_04Hr.adj.pvalue < 0.05 | Inf_08Hr.adj.pvalue < 0.05 | Inf_12Hr.adj.pvalue < 0.05 |  Inf_24Hr.adj.pvalue < 0.05 )
BP <- unique(Bouhaddou$Gene_Name)

# Annotate colors
net %v% "BP" = (levels %in% BP)

set.seed(1)
p_rb <-  ggnet2(net, label = BP, color = "BP",  legend.position = "none", mode = "fruchtermanreingold",
                layout.par = list(repulse.rad = 2000), alpha = 0.75, size = "size_of_node", size.cut = 10,
                palette = c("FALSE" = "grey", "TRUE" = "firebrick"),
                label.size = 2) +  guides( size = FALSE)
cowplot::ggsave2(p_rb, file = "../output/string/RAP-network_Bouhaddou.pdf", width = 4, height = 3)

