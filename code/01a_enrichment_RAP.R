library(data.table)
library(dplyr)
library(Matrix)
library(BuenColors)
library(network)

# Import database
if(!exists("string_v11")){
  source('00_string_v11.R')
}

# Import gene names
rap <- fread("../data/SCoV2_RAPms.txt") 

rap_filt <- rap%>%
  filter(adj.P.Val.SCoV2.over.RMRP < 0.2 & species == "HOMO SAPIENS" & logFC.SCoV2.over.RMRP > 0)
rap_genes <- rap_filt %>% pull(geneSymbol)

# Make string network
bulk <- fread("../data/SCoV2_bulkProteome.txt") %>% arrange(log(P.Value.SCoV.over.Mock)) %>%
  filter(species == "HOMO SAPIENS") 
all_ms_genes <- unique(bulk$geneSymbol, rap_genes)

observed <- pull_string_v11(min_score = 550, genes_keep = rap_genes) %>%
  filter(node1 != node2)

permuted_connections <- sapply(1:1000, function(i){
  set.seed(i)
  perm <- pull_string_v11(min_score = 550, genes_keep = sample(all_ms_genes, length(rap_genes))) %>%
    filter(node1 != node2)
  dim(perm)[1]
})

(dim(observed)[1]-mean(permuted_connections))/sd(permuted_connections)
mean(permuted_connections)

p1 <- ggplot(data.frame(connections = permuted_connections), aes(x = connections)) +
  geom_histogram(fill = "lightgrey", color = "black") + scale_x_log10() + 
  geom_vline(xintercept = 1534, color = "firebrick") +
  pretty_plot(fontsize = 8 ) + L_border() + scale_y_continuous(expand = c(0,0))+
  labs(x = "Total connections (log10 scale)", y = "Count")
cowplot::ggsave2(p1, file = "../output/permuted_RAP_string_network.pdf", width = 2.0, height = 1.8)
