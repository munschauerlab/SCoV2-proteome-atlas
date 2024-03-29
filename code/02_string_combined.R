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
all <- fread("../data/SCoV2_bulkProteome.txt") %>% arrange(log(P.Value.SCoV.over.Mock)) %>%
  filter(species == "HOMO SAPIENS") 

top1ms <- all[!duplicated(all$geneSymbol),]
rap <- fread("../data/SCoV2_RAPms.txt") %>%
  filter(adj.P.Val.SCoV2.over.RMRP < 0.05 & species == "HOMO SAPIENS" & logFC.SCoV2.over.RMRP > 0)

merge_df <- left_join(rap, top1ms, by = c("geneSymbol"))

top1ms <- top1ms %>% 
  filter(adj.P.Val.SCoV.over.Mock < 0.05)

bulk_down <- top1ms %>% filter(logFC.SCoV.over.Mock < 0) %>% pull(geneSymbol)

# Query
df <- pull_string_v11(min_score = 700, genes_keep = unique(c(rap$geneSymbol, top1ms$geneSymbol)))

# Reorientate
df_filt <- df[(df$node1 %in% top1ms$geneSymbol) | (df$node2 %in% top1ms$geneSymbol) ,]
df_filt$N1inRAP <- df_filt$node1 %in% rap$geneSymbol
df_filt$N2inRAP <- df_filt$node2 %in% rap$geneSymbol

df_filt2 <- df_filt[(df_filt$N1inRAP + df_filt$N2inRAP) ==1 ,]

write.xlsx(df_filt2, "../output/TableS7_interactome_regulated_network.xlsx",
           sheetName = "Network", 
           col.names = TRUE, row.names = FALSE, append = TRUE)

# Visualize network
edges <- df_filt2[,c("node1", "node2", "combined_score")]
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
  levels %in% c(rap$geneSymbol) & levels %in% c(top1ms$geneSymbol) ~ "both", 
  levels %in% c(rap$geneSymbol) ~ "RAP",
  
  levels %in% bulk_down ~ "BulkDown",
  TRUE ~ "BulkUp"
)
table(color_vec)
net %v% "color" = color_vec

i = 17
lapply(1:25, function(i){
  set.seed(i)
  pAI <- ggnet2(net, label = c(rap$geneSymbol), color = "color",  legend.position = "bottom", 
                alpha = 0.75, size = "degree",
                palette = c("BulkUp" = "lightgrey", "BulkDown"= "black", "both" = "red", "RAP" =  "dodgerblue", "highlight" = "green"),label.size = 2) +  
    guides( size = FALSE)
  cowplot::ggsave2(pAI, file = paste0("../output/string/String_combined_RAP+bulk_2July2020_v11_FDR05_seed",as.character(i),".pdf"), width = 8, height = 8)

})
