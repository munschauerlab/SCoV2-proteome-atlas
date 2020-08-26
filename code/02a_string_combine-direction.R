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

# Total bulk
table(unique(c(df_filt2$node1, df_filt2$node2)) %in% bulk_down)
table(top1ms$logFC.SCoV.over.Mock < 0)
