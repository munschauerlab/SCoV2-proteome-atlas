#' ---
#' title: "SARS CoV2 RAP-MS interaction network"
#' author: "Caleb Lareau, Mathias Munschauer"
#' date: "`r Sys.Date()`"
#' output: html_document
#' ---

#' Select a specifc node to visualize STRING interactions associated with the specific protein. 
#' 
#' 
#' Genes associated with a 'response to virus' pathway are shown in purple while other genes previous associated with viral RNA are shown in blue. Other significant genes at an adj. P < 0.2 are shown in grey from the RAP-MS experiment.
#+ cache = FALSE, message = FALSE, warning = FALSE, echo = FALSE, eval = TRUE,
# Import libraries

library(data.table)
library(dplyr)
library(Matrix)
library(networkD3)
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

# Visualize network
levels = sort(unique(c(edges[["node1"]], edges[["node2"]])))
edges_df <- data.frame(
  from = factor(edges$node1, levels = levels),
  to = factor(edges$node2, levels = levels),
  value = edges$combined_score / 1000)

# Make node annotation properties
n_size_vec = rap$logFC.SCoV2.over.RMRP; names(n_size_vec) <- as.character(rap$geneSymbol)
size_vec <- unname(n_size_vec[levels])
n_p_vec = rap$adj.P.Val.SCoV2.over.RMRP; names(n_p_vec) <- as.character(rap$geneSymbol)
p_vec <- unname(n_p_vec[levels])
rtv <- c("C19orf66", "DDX1", "LSM14A", "DDX3X", "PCBP2", "ACTA2", "CFL1")

viral_rna <- c("ATP1A1","CAPRIN1","CFL1","CSDE1","DDX1","DDX3X","EEF1A1","EEF2","EIF3E","EIF3G","EIF3L","EIF4B","EIF4G1","EIF4H","EIF5A","G3BP1","G3BP2",
               "HNRNPA1","HNRNPA2B1","IGF2BP1","IGF2BP2","LIN28B","LSM14A","MOV10","PABPC1","PCBP2","PEBP1","PFN1","PPIA","RPL13","RPL15","RPL18A","RPL21",
               "RPL28","RPL3","RPL36A","RPL6","RPL7A","RPL8","RPS11","RPS12","RPS14","RPS2","RPS26","RPS3","RPS4X","RPS5","SND1","SYNCRIP","UPF1","YBX1")

color_vec <- case_when(
  levels %in% rtv ~ "Response to Virus",
  levels %in% viral_rna ~ "Viral-associated",
  TRUE ~ "Other"
)
#table(color_vec)

nodes_df <- data.frame(
  id = levels,
  group = color_vec,
  value = size_vec * 10,
  title = paste0("<p>Node: ", levels, "<br>",
                 "log2FC: ",
                 as.character(round(size_vec, 2)),
                 "<br>",
                 "adj. P: ",
                 as.character(round(p_vec, 4)),
                 "</p>")
)

#+ cache = FALSE, message = FALSE, warning = FALSE, echo = FALSE, eval = TRUE, fig.width=10, fig.height=10
library(visNetwork)
visNetwork( nodes_df, edges_df) %>% 
  visIgraphLayout() %>% 
  visNodes(font = list(size = 40)) %>%
  visPhysics(stabilization = FALSE)  %>%
  visGroups(groupname = "Other", color = "grey") %>%
  visGroups(groupname = "Response to Virus", color = "#912cee") %>%
  visGroups(groupname = "Viral-associated", color = "#104e8b") %>%
  visEdges(color = list(highlight = "blue", hover = "blue",  smooth = FALSE)) %>%
  visOptions(highlightNearest = list(enabled = TRUE, degree = 1,
                                     labelOnly = FALSE, hover = TRUE),
             nodesIdSelection = TRUE) 
