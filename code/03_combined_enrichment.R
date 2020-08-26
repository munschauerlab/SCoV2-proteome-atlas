library(data.table)
library(dplyr)
library(Matrix)
library(BuenColors)

# Import database
if(!exists("string_v11")){
  source('00_string_v11.R')
}
# Import gene names
all <- fread("../data/SCoV2_bulkProteome.txt") %>% arrange(log(P.Value.SCoV.over.Mock)) %>%
  filter(species == "HOMO SAPIENS") 

top1ms <- all[!duplicated(all$geneSymbol),]


top1ms <- top1ms %>% 
  filter(adj.P.Val.SCoV.over.Mock < 0.05)

# Import RAP
rap <- fread("../data/SCoV2_RAPms.txt") %>%
  filter(adj.P.Val.SCoV2.over.RMRP < 0.2 & species == "HOMO SAPIENS" & logFC.SCoV2.over.RMRP > 0)
rap_genes <- rap %>% pull(geneSymbol)

String_df <- pull_string_v11(min_score = 700, genes_keep = unique(c(top1ms$geneSymbol)))



# Establish centrality
String_df_filt <- String_df[(String_df$node1 %in% top1ms$geneSymbol) & (String_df$node2 %in% top1ms$geneSymbol) ,]
count_df <- data.frame(
  gene = c(String_df_filt[[1]], String_df_filt[[2]])
) %>% group_by(gene) %>% summarize(count = n()) %>% mutate(inRAP= gene %in% rap_genes) %>%
  arrange(desc(count)) %>% mutate(rank = 1:n())


wilcox.test(count_df$count~count_df$inRAP)

count_df %>% group_by(inRAP) %>%
  summarize(mean(rank), mean(count), n = n())

