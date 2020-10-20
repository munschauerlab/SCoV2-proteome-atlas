library(data.table)
library(BuenColors)
library(annotables)
library(dplyr)
library(Matrix)
library(msigdbr)
library(tidyverse)

## import the MACS output
cnbp <- diffloop::bedToGRanges("../data/macs2/CNBP_ns_summits.bed") %>% diffloop::rmchr()
larp1 <- diffloop::bedToGRanges("../data/macs2/LARP1_ns_summits.bed")%>% diffloop::rmchr()
protein_coding_rna <- data.frame(grch38[,c(3:6,8)] %>% filter(biotype == "protein_coding")) %>%
  makeGRangesFromDataFrame(keep.extra.columns = TRUE, ignore.strand = TRUE)

cnbp_genes <- protein_coding_rna$symbol[unique(subjectHits(findOverlaps(cnbp, protein_coding_rna)))]
larp1_genes <- protein_coding_rna$symbol[unique(subjectHits(findOverlaps(larp1, protein_coding_rna)))]
all_genes <- protein_coding_rna$symbol

# Get msigdb
h_df_c5bp = msigdbr(species = "Homo sapiens") %>% dplyr::filter(gs_cat == "C5" & gs_subcat == "BP")
pathway_list_c5bp = h_df_c5bp %>% split(x = .$gene_symbol, f = .$gs_name)

computeFisherOverlapP <- function(clip_genes){
  
  # Do Fisher tests manually
  h_df_filt <- h_df_c5bp %>% filter(human_gene_symbol %in% all_genes)
  n_denom <- length(unique(h_df_filt$human_gene_symbol))
  
  # Function to get all of the 
  get_fisher_tests <- function(overlap_genes,  min_size = 3){
    h_df_filt2 <- h_df_filt %>% mutate(hit = human_gene_symbol %in% overlap_genes)
    analysis <- h_df_filt2 %>%group_by(gs_name) %>%
      summarize(a = sum(hit), n = n(), gene_hits =  paste(human_gene_symbol[hit],collapse = ",") ) %>%
      mutate(b = n-a, c = length(overlap_genes)-a) %>% mutate(d = n_denom - a - b -c) %>%
      filter(a >= min_size) %>%
      nest(-gs_name, -n, -gene_hits) %>% 
      mutate(matrix = map(data, ~matrix(unlist(.x), nrow = 2))) %>% 
      mutate(fisher = map(matrix, ~fisher.test(.x))) %>% 
      mutate(stats = map(fisher, ~broom::glance(.x)))
    
    if(dim(analysis)[1] == 0){
      return(NULL)
    }
    
    # unnest to report statistical significance
    analysis %>%
      unnest(stats) %>%
      dplyr::select(gs_name, gene_hits, n, p.value, odds = estimate) %>%
      mutate(padj = p.adjust(p.value)) %>% arrange(p.value)
    }
  get_fisher_tests(clip_genes)
}

# Top Filt uses the differential bulk proteome
# all uses all genes detected by the proteome
cnbp_enrichments <- computeFisherOverlapP(cnbp_genes)
larp1_enrichments <- computeFisherOverlapP(larp1_genes)

write.table(cnbp_enrichments, file = "../output/top_go_enrich-CNBP_clip.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

write.table(larp1_enrichments, file = "../output/top_go_enrich-LARP1_clip.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)

