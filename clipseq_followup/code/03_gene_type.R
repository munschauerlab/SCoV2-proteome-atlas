library(dplyr)
library(data.table)
library(BuenColors)
library(GenomicRanges)

## import the MACS output
keep_chrs <- paste0("chr", c(as.character(1:22), "X", "Y"))
cnbp <- diffloop::bedToGRanges("../data/macs2/CNBP_ns_summits.bed")
larp1 <- diffloop::bedToGRanges("../data/macs2/LARP1_ns_summits.bed")
cnbp <- cnbp[cnbp@seqnames %in% keep_chrs]
larp1 <- larp1[larp1@seqnames %in% keep_chrs]

# import RNA biotypes
rnatypes <- fread("../data/metagene_pileup/hg38_rnatypes.bed",
                  col.names = c("chr", "start", "end", "g", "m", "a", "what", "type"))
sort(table(rnatypes$type))
rt_gr <- makeGRangesFromDataFrame(data.frame(rnatypes), keep.extra.columns = TRUE)
protein_coding_rna <- rt_gr[rnatypes$type == "protein_coding"]
lncrna <- rt_gr[rnatypes$type == "lncRNA"]
pseudo <- rt_gr[grepl("pseudo", rnatypes$type)]

# Annotate further
x5utr <- diffloop::bedToGRanges("../data/metagene_pileup/hg38_merged_5utr.bed")
cds <- diffloop::bedToGRanges("../data/metagene_pileup/hg38_merged_cds.bed")
x3utr <- diffloop::bedToGRanges("../data/metagene_pileup/hg38_merged_3utr.bed")

# Function to score imported RNAs
score_assign_prioirity <- function(gr_in){
  protein_hits <- unique(queryHits(findOverlaps(gr_in,protein_coding_rna)))
  lncrna_hits <- unique(queryHits(findOverlaps(gr_in,lncrna)))
  pseudo_hits <- unique(queryHits(findOverlaps(gr_in,pseudo)))
  v <- 1:length(gr_in)
  out <- case_when(
    v %in% protein_hits ~ "zprotein_coding",
    v %in% lncrna_hits ~ "lncRNA",
    v %in% pseudo_hits ~ "pseudogene",
    TRUE~ "aother"
    
  )
}

larp1p <- table(score_assign_prioirity(larp1))/length(larp1) * 100
cnbpp <- table(score_assign_prioirity(cnbp))/length(cnbp) *100

p1 <- ggplot(data.frame(cnbpp), aes(x = 1, y = Freq, fill = Var1)) +
  geom_bar(stat = "identity", width = 0.2, color = "black", alpha = 0.5) +
  pretty_plot(fontsize=7) + L_border() +
  scale_x_continuous(limits = c(0.7, 1.3)) +
  scale_y_continuous(expand = c(0,0)) + labs(x = "CNBP", y = "% of CLIP peaks", fill = "") +
  theme(legend.position = "bottom") +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  scale_fill_manual(values = c("purple4", "firebrick", "darkgrey", "black"))
cowplot::ggsave2(p1, file = "../plots/stacked_bar_cnbp.pdf", width = 2.5, height = 2)

p2 <- ggplot(data.frame(larp1p), aes(x = 1, y = Freq, fill = Var1)) +
  geom_bar(stat = "identity", width = 0.2, color = "black", alpha = 0.5) +
  pretty_plot(fontsize=7) + L_border() +
  scale_x_continuous(limits = c(0.7, 1.3)) +
  scale_y_continuous(expand = c(0,0)) + labs(x = "LARP1", y = "% of CLIP peaks", fill = "") +
  theme(legend.position = "bottom") +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  scale_fill_manual(values = c("purple4", "firebrick", "darkgrey", "black"))
cowplot::ggsave2(p2, file = "../plots/stacked_bar_larp1.pdf", width = 2.5, height = 2)

score_assign_prioirity_within_coding <- function(gr_in){
  
  # Filter for protein overlapping
  protein_hits <- unique(queryHits(findOverlaps(gr_in,protein_coding_rna)))
  gr_in2 <- gr_in[1:length(gr_in) %in% protein_hits]
  
  ov5 <- unique(queryHits(findOverlaps(gr_in2,x5utr)))
  ovcds <- unique(queryHits(findOverlaps(gr_in2,cds)))
  ov3 <- unique(queryHits(findOverlaps(gr_in2,x3utr)))
  
  v <- 1:length(gr_in2)
  
  # Slight sdifferences in annotation complicate the overlap
  out <- case_when(
    v %in% ovcds ~ "CDS",
    v %in% ov5 ~ "5UTR",
    v %in% ov3 ~ "a3UTR",
    TRUE~ "none"
  )
  out[out!= "none"]
}

pct_vec <- function(v){
  table(v)/length(v) *100
}

cnbpp_z <- pct_vec(score_assign_prioirity_within_coding(cnbp)) 
larp1p_z <- pct_vec(score_assign_prioirity_within_coding(larp1))


p3 <- ggplot(data.frame(cnbpp_z), aes(x = 1, y = Freq, fill = v)) +
  geom_bar(stat = "identity", width = 0.2, color = "black", alpha = 0.5) +
  pretty_plot(fontsize=7) + L_border() +
  scale_x_continuous(limits = c(0.7, 1.3)) +
  scale_y_continuous(expand = c(0,0)) + labs(x = "CNBP", y = "% of CLIP coding peaks", fill = "") +
  theme(legend.position = "bottom") +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  scale_fill_manual(values = c("green4", "darkblue","darkgrey"))
cowplot::ggsave2(p3, file = "../plots/stacked_bar_cnbp_codeanno.pdf", width = 2.5, height = 2)

p4 <- ggplot(data.frame(larp1p_z), aes(x = 1, y = Freq, fill = v)) +
  geom_bar(stat = "identity", width = 0.2, color = "black", alpha = 0.5) +
  pretty_plot(fontsize=7) + L_border() +
  scale_x_continuous(limits = c(0.7, 1.3)) +
  scale_y_continuous(expand = c(0,0)) + labs(x = "LARP1", y = "% of CLIP coding peaks", fill = "") +
  theme(legend.position = "bottom") +
  theme(axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  scale_fill_manual(values = c("green4",  "darkblue", "darkgrey"))
cowplot::ggsave2(p4, file = "../plots/stacked_bar_larp1_codeanno.pdf", width = 2.5, height = 2)


