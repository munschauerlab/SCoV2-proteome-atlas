library(data.table)
library(dplyr)

# Import database
if(!exists("string_v11")){
  source('00_string_v11.R')
}
# Import gene names
all <- fread("../data/all_proteomics.txt") %>% arrange(desc(Log.P.Value.SCoV.over.Mock)) %>%
  filter(adj.P.Val.SCoV.over.Mock < 0.05 & geneSymbol %ni% viral_genes)
top1ms <- all[!duplicated(all$geneSymbol),]
rap <- fread("../data/MM_proteomics_data_limmaFDR.txt") %>%
  filter(adj.P.Val.SARS.CoV.over.RMRP < 0.2 & geneSymbol %ni% viral_genes)

# Query
df <- pull_string_v11(min_score = 700, genes_keep = unique(c(rap$geneSymbol, all$geneSymbol)))

# Reorientate
df_filt <- df[(df$node1 %in% top1ms$geneSymbol) | (df$node2 %in% top1ms$geneSymbol) ,]
df_filt$N1inRAP <- df_filt$node1 %in% top1ms$geneSymbol
df_filt$N2inRAP <- df_filt$node2 %in% top1ms$geneSymbol

df_filt2 <- df_filt[df_filt$N1inRAP + df_filt$N2inRAP < 2,]

#write.table(df_filt2, file = "../output/string/RAP_MS_STRING_network.tsv", 
#            sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

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
  levels %in% c(rap$geneSymbol) & levels %in% c(all$geneSymbol) ~ "both", 
  levels %in% c(rap$geneSymbol) ~ "RAP",
  TRUE ~ "Bulk"
)
net %v% "color" = color_vec

set.seed(3)
pAI <- ggnet2(net, label = rap$geneSymbol, color = "color",  legend.position = "bottom", 
              palette = c("Bulk" = "grey", "both" = "red", "RAP" =  "dodgerblue"), size = 3, label.size = 4)
cowplot::ggsave2(pAI, file = "../output/string/String_combined_RAP+bulk_27June2020_v11.pdf", width = 8, height = 8)


