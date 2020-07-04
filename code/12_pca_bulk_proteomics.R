library(BuenColors)
library(dplyr)
library(data.table)

# Import differential gene expression
all <- data.frame(fread("../data/bulk_proteome_raw_data.txt"))
data.frame(prcomp(t(data.matrix(all[,c(1:6)])), scale = TRUE, center = TRUE)$x,
           what = c("Mock", "Mock", "Mock", "SCoV", "SCoV", "SCoV")) -> pldf

ggplot(pldf, aes(x = PC1, y = PC2, color = what)) +
  geom_point() + pretty_plot(fontsize = 10) +
  L_border() + scale_color_manual(values = c("dodgerblue3", "firebrick")) +
  labs(color = "") + theme(legend.position = "bottom") -> p1

cowplot::ggsave2(p1, file = "../output/PCAplot_bulk.pdf", width = 2, height = 2.4)
