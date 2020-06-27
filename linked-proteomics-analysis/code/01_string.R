library(STRINGdb)
library(data.table)
library(dplyr)

# Import database
string_db <- STRINGdb$new( version="10", species=9606, score_threshold=250, input_directory="")

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

write.table(df_filt2, file = "../RAP_MS_STRING_network.tsv", 
            sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)

dim(df_filt)
dim(df_filt2)
