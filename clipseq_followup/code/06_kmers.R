library(rtracklayer)         
library(BSgenome)                  
library(BSgenome.Hsapiens.UCSC.hg38)
library(Biostrings)                   
library(dplyr)
library(motifmatchr)
library(chromVAR)
library(universalmotif)
library(SummarizedExperiment)
library(BuenColors)

cnbp <-  import("../data/macs2/CNBP_ns_peaks.narrowPeak") %>% reduce()
cnbp <- cnbp[cnbp@seqnames!= "NC_045512.2"]

larp1 <-  import("../data/macs2/LARP1_ns_peaks.narrowPeak")%>% reduce()
larp1 <- larp1[larp1@seqnames!= "NC_045512.2"]
x5utr <- diffloop::bedToGRanges("../data/metagene_pileup/hg38_merged_5utr.bed")
larp1 <- larp1[unique(queryHits(findOverlaps(larp1, x5utr)))]

cnbp_km <- matchKmers(5, cnbp, BSgenome.Hsapiens.UCSC.hg38)
cnbp_permute <- matchKmers(5, shuffle_sequences(getSeq(BSgenome.Hsapiens.UCSC.hg38,cnbp)))

cnbp_enrich <- colSums(assays(cnbp_km)[["matches"]])/colSums(assays(cnbp_permute)[["matches"]])
tail(sort(cnbp_enrich))

larp1_km <- matchKmers(5, larp1, BSgenome.Hsapiens.UCSC.hg38)
larp1_permute <- matchKmers(5, shuffle_sequences(getSeq(BSgenome.Hsapiens.UCSC.hg38,larp1)))

larp1_enrich <- colSums(assays(larp1_km)[["matches"]])/colSums(assays(larp1_permute)[["matches"]])

larp1_df <- data.frame(
  rank = rev(1:length(larp1_enrich)),
  value = unname(sort(larp1_enrich)),
  name = names(sort(larp1_enrich))
) %>% arrange(rank)
larp1_df$rc <- as.character(reverseComplement(DNAStringSet(larp1_df$name)))
head(larp1_df, 20)


cnbp_df <- data.frame(
  rank = rev(1:length(cnbp_enrich)),
  value = unname(sort(cnbp_enrich)),
  name = names(sort(cnbp_enrich))
) %>% arrange(rank)
cnbp_df$rc <- as.character(reverseComplement(DNAStringSet(cnbp_df$name)))

head(cnbp_df, 10)
