library(data.table)
library(dplyr)
library(Matrix)
library(BuenColors)


# Import gene names
all <- fread("../data/SCoV2_bulkProteome.txt") %>% arrange(log(P.Value.SCoV.over.Mock)) %>%
  filter(species == "HOMO SAPIENS") 

rap <- fread("../data/SCoV2_RAPms.txt") %>%
  filter(species == "HOMO SAPIENS")

mdf <- merge(rap, all, by = "id")

cor(mdf$logFC.SCoV2.over.RMRP, mdf$logFC.SCoV.over.Mock)
p1 <- ggplot(mdf, aes(x = logFC.SCoV2.over.RMRP, y = logFC.SCoV.over.Mock)) +
  geom_point(size = 0.5) + pretty_plot(fontsize = 8) + L_border() +
  geom_smooth(method='lm', formula= y~x) + 
  labs(x = "logFC RAP-MS Fold Change", y = "logFC Bulk Fold Change")
cowplot::ggsave2(p1, file = "../output/scatter_RAP_bulk.pdf", width = 2, height = 2)

summary(lm(mdf$logFC.SCoV2.over.RMRP~ mdf$logFC.SCoV.over.Mock))
