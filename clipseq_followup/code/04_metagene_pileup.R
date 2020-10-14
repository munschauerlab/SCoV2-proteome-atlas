library(data.table)
library(ggplot2)
library(BuenColors)

CNBP_cds <- fread("../data/metagene_pileup/CNBP_cds.tsv.gz", header = FALSE)
CNBP_5utr <- fread("./data/metagene_pileup/CNBP_5utr.tsv.gz", header = FALSE)
CNBP_3utr <- fread("./data/metagene_pileup/CNBP_3utr.tsv.gz", header = FALSE)

LARP1_cds <- fread("./data/metagene_pileup/LARP1_cds.tsv.gz", header = FALSE)
LARP1_5utr <- fread("./data/metagene_pileup/LARP1_5utr.tsv.gz", header = FALSE)
LARP1_3utr <- fread("./data/metagene_pileup/LARP1_3utr.tsv.gz", header = FALSE)

larp1_df <- rbind(
  data.frame(
    position = 1:100, 
    score = colMeans(data.matrix(data.frame(LARP1_5utr[,c(7:106)])), na.rm = TRUE),
    what = "5UTR"
  ),
  data.frame(
    position = 101:600, 
    score = colMeans(data.matrix(data.frame(LARP1_cds[,c(7:506)])), na.rm = TRUE),
    what = "CDS"
  ),
  data.frame(
    position = 601:700, 
    score = colMeans(data.matrix(data.frame(LARP1_3utr[,c(7:106)])), na.rm = TRUE),
    what = "3UTR"
  ))


cnbp_df <- rbind(
  data.frame(
    position = 1:100, 
    score = colMeans(data.matrix(data.frame(CNBP_5utr[,c(7:106)])), na.rm = TRUE),
    what = "5UTR"
  ),
  data.frame(
    position = 101:600, 
    score = colMeans(data.matrix(data.frame(CNBP_cds[,c(7:506)])), na.rm = TRUE),
    what = "CDS"
  ),
  data.frame(
    position = 601:700, 
    score = colMeans(data.matrix(data.frame(CNBP_3utr[,c(7:106)])), na.rm = TRUE),
    what = "3UTR"
  ))

p1 <- ggplot(cnbp_df, aes(x = position, y = score, color = what)) +
  geom_vline(xintercept = c(100, 600), linetype = 2) +
  geom_line() +
  pretty_plot() + L_border() +
  theme(legend.position = "bottom") +
  scale_color_manual(values = c("firebrick", "purple3", "dodgerblue4")) +
  labs(color = "", y = "CNBP signal") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

p2 <- ggplot(larp1_df, aes(x = position, y = score, color = what)) +
  geom_vline(xintercept = c(100, 600), linetype = 2) +
  geom_line() +
  pretty_plot() + L_border() +
  theme(legend.position = "bottom") +
  scale_color_manual(values = c("firebrick", "purple3", "dodgerblue4")) +
  labs(color = "", y = "LARP1 signal") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

cowplot::ggsave2(p1, file = "../plots/CNBP_metagene.pdf", width = 3, height = 2)
cowplot::ggsave2(p2, file = "../plots/LARP1_metagene.pdf", width = 3, height = 2)
