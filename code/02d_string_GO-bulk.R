library(data.table)
library(dplyr)
library(Matrix)
library(network)
library(msigdbr)
library(tidyverse)

# Get msigdb
h_df_c5bp = msigdbr(species = "Homo sapiens") %>% dplyr::filter(gs_cat == "C5" & gs_subcat == "BP")
pathway_list_c5bp = h_df_c5bp %>% split(x = .$gene_symbol, f = .$gs_name)

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

top1ms_filt <- top1ms %>% 
  filter(adj.P.Val.SCoV.over.Mock < 0.05)

bulk_down <- top1ms_filt %>% filter(logFC.SCoV.over.Mock < 0) %>% pull(geneSymbol)

# Query
df <- pull_string_v11(min_score = 700, genes_keep = unique(c(rap$geneSymbol, top1ms_filt$geneSymbol)))

# Reorientate
df_filt <- df[(df$node1 %in% top1ms$geneSymbol) | (df$node2 %in% top1ms_filt$geneSymbol) ,]
df_filt$N1inRAP <- df_filt$node1 %in% rap$geneSymbol
df_filt$N2inRAP <- df_filt$node2 %in% rap$geneSymbol

df_filt2 <- df_filt[(df_filt$N1inRAP + df_filt$N2inRAP) ==1 ,]

# Do Fisher tests manually
h_df_filt <- h_df_c5bp %>% filter(human_gene_symbol %in% top1ms_filt$geneSymbol)
n_denom <- length(unique(h_df_filt$human_gene_symbol))

# Function to get all of the 
get_fisher_tests <- function(overlap_genes, seed_gene = "", min_size = 3){
  h_df_filt2 <- h_df_filt %>% mutate(hit = human_gene_symbol %in% overlap_genes)
  analysis <- h_df_filt2 %>%group_by(gs_name) %>%
    summarize(a = sum(hit), n = n()) %>%
    mutate(b = n-a, c = length(overlap_genes)-a) %>% mutate(d = n_denom - a - b -c) %>%
    filter(a >= min_size) %>%
    nest(-gs_name, -n) %>% 
    mutate(matrix = map(data, ~matrix(unlist(.x), nrow = 2))) %>% 
    mutate(fisher = map(matrix, ~fisher.test(.x))) %>% 
    mutate(stats = map(fisher, ~broom::glance(.x)))
  
  if(dim(analysis)[1] == 0){
    return(NULL)
  }
  
  # unnest to report statistical significance
  analysis %>%
    unnest(stats) %>%
    dplyr::select(gs_name, n, p.value, odds = estimate) %>%
    mutate(padj = p.adjust(p.value)) %>% arrange(p.value) %>%
    filter(padj < 0.01) %>% mutate(seed_gene = seed_gene)
}
rap_genes <- rap$geneSymbol

lapply(rap_genes, function(gene){
  print(gene)
  seed_gene_connections <- c(df_filt2[df_filt2$node1==gene, ][["node2"]],df_filt2[df_filt2$node2==gene, ][["node1"]])
  get_fisher_tests(seed_gene_connections, gene)
}) %>% rbindlist() %>% data.frame() -> all_community_enrichments

top_community_enrichments <- all_community_enrichments %>% group_by(seed_gene) %>% arrange(seed_gene,padj) %>% top_n(5, wt =-log10(padj))
write.table(top_community_enrichments, file = "../output/top_go_enrich_network_communities.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
