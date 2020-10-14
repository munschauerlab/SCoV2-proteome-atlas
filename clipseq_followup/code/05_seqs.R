library(rtracklayer)                  # for import()
library(BSgenome)                     # for getSeq()
library(BSgenome.Hsapiens.UCSC.hg38)
library(Biostrings)                   # for writeXStringSet()

cnbp <-  import("../data/macs2/CNBP_ns_peaks.narrowPeak")
cnbp_seq <- getSeq(BSgenome.Hsapiens.UCSC.hg38,cnbp[cnbp@seqnames!= "NC_045512.2"])

larp1 <-  import("../data/macs2/LARP1_ns_peaks.narrowPeak")
larp1_seq <- getSeq(BSgenome.Hsapiens.UCSC.hg38,larp1[larp1@seqnames!= "NC_045512.2"])

names(cnbp_seq) <- paste0("c", as.character(1:length(cnbp_seq)))
names(larp1_seq) <- paste0("l", as.character(1:length(larp1_seq)))
writeXStringSet(cnbp_seq, "../output/seqs_cnbp.fasta")
writeXStringSet(larp1_seq, "../output/seqs_larp1.fasta")
