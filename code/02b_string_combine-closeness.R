library(data.table)
library(dplyr)
library(Matrix)
library(igraph)
library(network)
library(BuenColors)
"%ni%" <- Negate("%in%")

# Import database
if(!exists("string_v11")){
  source('00_string_v11.R')
}
# Import gene names
all <- fread("../data/SCoV2_bulkProteome.txt") %>% arrange(log(P.Value.SCoV.over.Mock)) %>%
  filter(species == "HOMO SAPIENS") 

top1ms <- all[!duplicated(all$geneSymbol),]
rap <- fread("../data/SCoV2_RAPms.txt") %>%
  filter(adj.P.Val.SCoV2.over.RMRP < 0.2 & species == "HOMO SAPIENS" & logFC.SCoV2.over.RMRP > 0)

merge_df <- left_join(rap, top1ms, by = c("geneSymbol"))

top1ms_filt <- top1ms %>% 
  filter(adj.P.Val.SCoV.over.Mock < 0.05)

# Query
df <- pull_string_v11(min_score = 700, genes_keep = unique(c(rap$geneSymbol, top1ms$geneSymbol)))
levels = sort(unique(c(df$node1, df$node2)))
factor_df <- rbind(data.frame(
  n1 = factor(as.character(df$node1), levels = levels),
  n2 = factor(as.character(df$node2), levels = levels)
), data.frame(n1 = levels, n2 = levels))
igraph <- graph_from_edgelist(data.matrix(factor_df), directed = FALSE)
s.path<- shortest.paths(igraph, algorithm = "dijkstra")

# Make some summary statistics
rownames(s.path) <- levels
colnames(s.path) <- levels
mdf <- reshape2::melt(s.path)
mdf$in_rap1 <- mdf$Var1 %in% rap$geneSymbol
mdf$in_rap2 <- mdf$Var2 %in% rap$geneSymbol
mdf$bulk1 <- mdf$Var1 %in% top1ms_filt$geneSymbol
mdf$bulk2 <- mdf$Var2 %in% top1ms_filt$geneSymbol
mdf_filt <- mdf %>% filter((in_rap1 | in_rap2) & (bulk1 | bulk2))

min_df <- data.frame(
  gene = levels[as.numeric(c(mdf_filt$Var1,mdf_filt$Var2))],
  value = c(mdf_filt$value,mdf_filt$value)
) %>% group_by(gene) %>% summarize(connection = min(value)) 

p1 <- min_df %>% group_by(connection) %>% summarize(count =n()) %>%
  mutate(prop = count/ sum(count)*100) %>%
  ggplot(aes(x = as.factor(as.character(connection)), y = prop)) + 
  geom_histogram(stat = "identity", color = "black", fill = "lightgrey") +
  pretty_plot(fontsize = 8) + L_border() +
  scale_y_continuous(expand = c(0,0)) + labs(x = "Degrees of sep.", y = "% of differential proteins")
cowplot::ggsave2(p1, file = "../output/histogram_connectivity.pdf", 
                width = 1.8, height = 1.8)
