library(STRINGdb)
library(data.table)
library(dplyr)

# Import database
string_db <- STRINGdb$new( version="10", species=9606, score_threshold=250, input_directory="../string/")

# Import gene names
all <- fread("../data/all_proteomics.txt") %>% arrange(desc(Log.P.Value.SCoV.over.Mock)) %>%
  filter(adj.P.Val.SCoV.over.Mock < 0.05)
top1ms <- all[!duplicated(all$geneSymbol),]
rap <- fread("../data/MM_proteomics_data_limmaFDR.txt") %>%
  filter(adj.P.Val.SARS.CoV.over.RMRP < 0.2)

# Map to STRING Ids
top1ms_mapped <- string_db$map( top1ms, "geneSymbol", removeUnmappedRows = TRUE )
rap_mapped <- string_db$map( rap, "geneSymbol", removeUnmappedRows = TRUE ) 


# Build interaction network
query_proteins <- c(top1ms_mapped$STRING_id, rap_mapped$STRING_id)
trans_vec <- c(top1ms_mapped$geneSymbol, rap_mapped$geneSymbol); names(trans_vec) <- query_proteins
df <- string_db$get_interactions(unique(query_proteins))[,c("from", "to", "combined_score")] %>%
  arrange(desc(combined_score))
colnames(df) <- c("A", "B", "score")
df$A_gene <- trans_vec[as.character(df$A)]
df$B_gene <- trans_vec[as.character(df$B)]

# Reorientate
df_filt <- df[(df$A %in% top1ms_mapped$STRING_id) | (df$B %in% top1ms_mapped$STRING_id) ,]
df_filt$AinRAP <- df_filt$A %in% top1ms_mapped$STRING_id
df_filt$BinRAP <- df_filt$B %in% top1ms_mapped$STRING_id

df_filt2 <- df_filt[df_filt$AinRAP + df_filt$BinRAP < 2,]

#write.table(df_filt2, file = "../output/string/RAP_MS_STRING_network.tsv", 
#            sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

# Visualize network
edges <- df_filt2[,c("A_gene", "B_gene", "score")]
colnames(edges) <- c("node1", "node2", "score")
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
  levels %in% c(rap$geneSymbol) & levels %in% c(all$geneSymbol) ~ "both", 
  levels %in% c(rap$geneSymbol) ~ "RAP",
  TRUE ~ "Bulk"
)
net %v% "color" = color_vec

set.seed(3)
pAI <- ggnet2(net, label = rap$geneSymbol, color = "color",  legend.position = "bottom", 
              palette = c("Bulk" = "grey", "both" = "red", "RAP" =  "dodgerblue"), size = 3, label.size = 4)
cowplot::ggsave2(pAI, file = "../output/string/String_combined_RAP+bulk_27June2020.pdf", width = 8, height = 8)


