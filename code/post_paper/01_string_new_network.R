library(data.table)
library(dplyr)
library(Matrix)
library(GGally)
library(network)

# Import database
if(!exists("string_v11")){
  source('00_string_v11_new.R')
}

# Import gene names
rap_genes <- fread("mm_new_genes.txt", header = FALSE) %>% pull(V1)

# Make string network
edges <- pull_string_v11(min_score = 800, genes_keep = rap_genes)
#write.xlsx(edges, "../output/TableS4_interactome_only_network.xlsx",
#           sheetName = "Interactome", 
#           col.names = TRUE, row.names = FALSE, append = TRUE)

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
all <- readxl::read_xlsx("Three_screen_data_long.xlsx")
top <- all %>% group_by(geneSymbol) %>% summarize(value = max(logFC.All))
n_size_vec = top$value; names(n_size_vec) <- as.character(top$geneSymbol)
size_vec <- unname(n_size_vec[levels])
net %v% "size_of_node" = size_vec

set.seed(1)

pAI <- ggnet2(net, label = levels, color = "color",  legend.position = "none", mode = "fruchtermanreingold", # [color_vec != "Other"]
              layout.par = list(repulse.rad = 2000), alpha = 0.75, size = "size_of_node", size.cut = 10,
              palette = c("Other" = "grey", "RP" =  "dodgerblue4", "viral_rna_bind" = "purple2"),
              label.size = 2) +  guides( size = FALSE)
cowplot::ggsave2(pAI, file = paste0("plots/String_RAPnetwork.13May2021.pdf"), width = 4, height = 3.5)



drugs <- fread("../../data/RAP-MS_dgidb_export_2020-06-29.tsv")
druggable <- drugs %>% pull(search_term) %>% unique()

net %v% "drug_target" = levels %in% druggable

set.seed(1)
pDrug <- ggnet2(net, label = druggable, color = "drug_target",  legend.position = "none", mode = "fruchtermanreingold",
              layout.par = list(repulse.rad = 2000), alpha = 0.75, size = "size_of_node", size.cut = 10,
              palette = c("FALSE" = "grey", "TRUE" = "firebrick"),
              label.size = 2) +  guides( size = FALSE)
cowplot::ggsave2(pDrug, file = paste0("plots/String_RAPnetwork_drugTargets.pdf"), width = 4, height = 3.5)


rb <- c("RPS2", "RPS14", "RPL15", "RPS10","NUDT3", "RPL6", "RPS3", "RPL13", "RPS4X", "RPS26", "RPS11", "RPL21", "RPS12", "RPS5", "RPL8", "RPL18A", "RPL3", "RPL36A", "RPL28", "RPL7A")
rtv <- c("C19orf66", "DDX1", "LSM14A", "DDX3X", "PCBP2", "ACTA2", "CFL1")
vp <- c("RAB6A", "RAB7B", "RAB1A", "EGFR", "PPIA", "SND1", "APOE", "HNRNPA1", "PCBP2", "DDX3X", "SYNCRIP", "EIF4G1", "EIF3G", "EIF4H", "EIF3H")
iep <- c("RAB6A", "RAB7B", "STOM", "C19orf66", "LSM14A", "PCBP2", "DDX3X", "ACTR2", "ACTB", "PPIA", "EEF2", "EEF1A1", "GPI")
vs <- vesicle

# Annotate colors
net %v% "rb" = (levels %in% rb)
net %v% "rtv" = (levels %in% rtv)
net %v% "vp" = (levels %in% vp)
net %v% "iep" = (levels %in% iep)
net %v% "vs" = (levels %in% vs)


set.seed(1)
p_rb <-  ggnet2(net, label = rb, color = "rb",  legend.position = "none", mode = "fruchtermanreingold",
                layout.par = list(repulse.rad = 2000), alpha = 0.75, size = "size_of_node", size.cut = 10,
                palette = c("FALSE" = "grey", "TRUE" = "firebrick"),
                label.size = 2) +  guides( size = FALSE)
cowplot::ggsave2(p_rb, file = "plots/network_rb.pdf", width = 4, height = 3)

set.seed(1)
p_vs <-  ggnet2(net, label = vs, color = "vs",  legend.position = "none", mode = "fruchtermanreingold",
                layout.par = list(repulse.rad = 2000), alpha = 0.75, size = "size_of_node", size.cut = 10,
                palette = c("FALSE" = "grey", "TRUE" = "firebrick"),
                label.size = 2) +  guides( size = FALSE)
cowplot::ggsave2(p_vs, file = "plots/network_vs.pdf", width = 4, height = 3)


set.seed(1)
p_rtv <-  ggnet2(net, label = rtv, color = "rtv",  legend.position = "none", mode = "fruchtermanreingold",
                 layout.par = list(repulse.rad = 2000), alpha = 0.75, size = "size_of_node", size.cut = 10,
                 palette = c("FALSE" = "grey", "TRUE" = "firebrick"),
                 label.size = 2) +  guides( size = FALSE)
cowplot::ggsave2(p_rtv, file = "plots/network_rtv.pdf", width = 4, height = 3)

set.seed(1)
p_vp <-  ggnet2(net, label = vp, color = "vp",  legend.position = "none", mode = "fruchtermanreingold",
                layout.par = list(repulse.rad = 2000), alpha = 0.75, size = "size_of_node", size.cut = 10,
                palette = c("FALSE" = "grey", "TRUE" = "firebrick"),
                label.size = 2) +  guides( size = FALSE)
cowplot::ggsave2(p_vp, file = "plots/network_vb.pdf", width = 4, height = 3)

set.seed(1)
p_iep <-  ggnet2(net, label = iep, color = "iep",  legend.position = "none", mode = "fruchtermanreingold",
                 layout.par = list(repulse.rad = 2000), alpha = 0.75, size = "size_of_node", size.cut = 10,
                 palette = c("FALSE" = "grey", "TRUE" = "firebrick"),
                 label.size = 2) +  guides( size = FALSE)
cowplot::ggsave2(p_iep, file = "plots/network_iep.pdf", width = 4, height = 3)



# Make new
Bouhaddou <- readxl::read_xlsx("../../data/Bouhaddou2020_TS1.xlsx") %>%
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

library(cowplot)
cowplot::ggsave2(plot_grid(pAI,plot_grid(p_rb, pDrug, ncol = 1), nrow = 1, rel_widths = c(1.25, 0.75)),
                 file = "plots/networks_3_grid_rev.pdf", width = 8, height = 5)


