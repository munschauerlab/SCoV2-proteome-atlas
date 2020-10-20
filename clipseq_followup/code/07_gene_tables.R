library(data.table)
library(BuenColors)
library(annotables)
library(dplyr)
library(GenomicRanges)
library(rtracklayer)
library(tidyverse)
library(xlsx)

## import the MACS output
cnbp <-  diffloop::bedToGRanges("../data/macs2/CNBP_ns_peaks.narrowPeak") %>% GenomicRanges::reduce() %>% diffloop::rmchr()
cnbp <- cnbp[cnbp@seqnames!= "NC_045512.2"]

larp1 <-  diffloop::bedToGRanges("../data/macs2/LARP1_ns_peaks.narrowPeak")%>% GenomicRanges::reduce()%>% diffloop::rmchr()
larp1 <- larp1[larp1@seqnames!= "NC_045512.2"]

# Annotate further
x5utr <- diffloop::bedToGRanges("../data/metagene_pileup/hg38_merged_5utr.bed") %>% diffloop::rmchr()
cds <- diffloop::bedToGRanges("../data/metagene_pileup/hg38_merged_cds.bed")%>% diffloop::rmchr()
x3utr <- diffloop::bedToGRanges("../data/metagene_pileup/hg38_merged_3utr.bed") %>% diffloop::rmchr()

protein_coding_rna <- data.frame(grch38[,c(3:6,8)] %>% filter(biotype == "protein_coding")) %>%
  makeGRangesFromDataFrame(keep.extra.columns = TRUE, ignore.strand = TRUE)

assign_within_coding <- function(gr_in){
  
  # Filter for protein overlapping
  ov <- findOverlaps(gr_in,protein_coding_rna)
  v <- queryHits(ov)
  d3 <- data.frame(diffloop::addchr(gr_in[v]))[,c(1:3)]
  
  df <- data.frame(
    peak = paste0(d3[[1]], ":", d3[[2]], "-", d3[[3]]),
    
    gene = protein_coding_rna$symbol[subjectHits(ov)]
  ) %>% unique()
  dfa <- aggregate(peak ~ gene, df, paste, collapse = "|")
  
  return(dfa)
}

larp1_df <- assign_within_coding(larp1)
cnbp_df <- assign_within_coding(cnbp)

# Annotate with published results
cnbp_pub <- readxl::read_xlsx("../data/pub/Benhalevy_CNBP_CLIP_genes.xlsx")[[13]]
mTOR_LARP1 <- readxl::read_xlsx("../data/pub/mTOR_LARP1 regulated_pnas.1912864117.sd03-1.xlsx", sheet = 3)[[1]]
cnbp_df$Benhalevy_CNBP_genes <- cnbp_df$gene %in% cnbp_pub
larp1_df$mTOR_LARP1_genes <- larp1_df$gene %in% mTOR_LARP1

table(cnbp_df$Benhalevy_CNBP_genes)
table(larp1_df$mTOR_LARP1_genes )

write.xlsx(cnbp_df, "../output/TableSX_CNBP_clip_genes.xlsx",
           sheetName = "CNBP_clip", 
           col.names = TRUE, row.names = FALSE, append = TRUE)

write.xlsx(larp1_df, "../output/TableSX_LARP1_clip_genes.xlsx",
           sheetName = "LARP1_clip", 
           col.names = TRUE, row.names = FALSE, append = TRUE)
