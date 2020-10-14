library(ChIPpeakAnno)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(BuenColors)

## import the MACS output
cnbp <- diffloop::bedToGRanges("CNBPmacs2_summits.bed")
larp1 <- diffloop::bedToGRanges("LARP1macs2_summits.bed")
ordering <- c("fiveUTRs", "threeUTRs", "Promoters", "Exons", "Introns",
  "immediateDownstream" )
aCNBP<-assignChromosomeRegion(cnbp, nucleotideLevel=FALSE, 
                              precedence=ordering, 
                              TxDb=TxDb.Hsapiens.UCSC.hg38.knownGene)
aLARP1<-assignChromosomeRegion(larp1, nucleotideLevel=FALSE, 
                               precedence=ordering, 
                               TxDb=TxDb.Hsapiens.UCSC.hg38.knownGene)

df <- rbind(data.frame(Freq=aCNBP$percentage, factor = "CNBP", what = names(aCNBP$percentage)), 
            data.frame(Freq=aLARP1$percentage, factor = "LARP1", what = names(aCNBP$percentage)))


orderingp <- c("fiveUTRs",  "Promoters", "Exons", "Introns","threeUTRs",
              "immediateDownstream",  "Intergenic.Region")
df$what_order <- factor(as.character(df$what), orderingp)

p1 <- ggplot(df,aes(x = what_order, y = Freq.Freq, fill = factor)) + 
  geom_bar(stat = "identity", position = 'dodge', color = 'black') +
  pretty_plot(fontsize = 7) + L_border() + theme(legend.position = "bottom") +
  scale_y_continuous(expand  = c(0,0)) +
  scale_fill_manual(values = c("lightgrey", "darkgrey")) +
  labs(x = "", y = "% CLIP peaks", color = "")
cowplot::ggsave2(p1, file = "overlaps_11October2020.pdf", width = 2.5, height = 2.7)
