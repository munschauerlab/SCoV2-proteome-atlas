library(rtracklayer)                  # for import()
library(BSgenome)                     # for getSeq()
library(BSgenome.Hsapiens.UCSC.hg38)
library(Biostrings)                   # for writeXStringSet()
library(dplyr)
"%ni%" <- Negate("%in%")

# Pull Larp1
larp1 <-  import("../data/macs2/LARP1_ns_peaks.narrowPeak")%>% reduce()
larp1 <- larp1[larp1@seqnames!= "NC_045512.2"]
x5utr <- diffloop::bedToGRanges("../data/metagene_pileup/hg38_merged_5utr.bed")
x5utr_neg <-x5utr[x5utr@strand == "-"]
x5utr_pos <-x5utr[x5utr@strand == "+"]

# Find stranded overlaps
pos_hits <- unique(queryHits(findOverlaps(larp1, x5utr_pos)))
neg_hits <- unique(queryHits(findOverlaps(larp1, x5utr_neg)))
both <- intersect(neg_hits, pos_hits)

larp1_pos <- larp1[pos_hits[pos_hits %ni% both]]
larp1_neg <- larp1[neg_hits[neg_hits %ni% both]]
larp1_pos@strand <- Rle("+", lengths = length(larp1_pos))
larp1_neg@strand <- Rle("-", lengths = length(larp1_neg))

larp1_seq <- c(getSeq(BSgenome.Hsapiens.UCSC.hg38, larp1_pos),
               getSeq(BSgenome.Hsapiens.UCSC.hg38, larp1_neg))
names(larp1_seq) <- paste0("l", as.character(1:length(larp1_seq)))
writeXStringSet(larp1_seq, "../output/seqs_larp1_stranded_5utrs.fasta")

# CNBP
cnbp <-  import("../data/macs2/CNBP_ns_peaks.narrowPeak") %>% reduce()
cnbp_seq <- getSeq(BSgenome.Hsapiens.UCSC.hg38,cnbp[cnbp@seqnames!= "NC_045512.2"])
names(cnbp_seq) <- paste0("c", as.character(1:length(cnbp_seq)))
writeXStringSet(cnbp_seq, "../output/seqs_cnbp.fasta")
