library(data.table)
library(dplyr)
library(Matrix)
library(BuenColors)

# Import gene names
all <- fread("../data/SCoV2_bulkProteome.txt") %>% arrange(log(P.Value.SCoV.over.Mock)) %>%
  filter(species == "HOMO SAPIENS") 

top1ms <- all[!duplicated(all$geneSymbol),]
top1schmidt <- top1ms[,c("geneSymbol", "logFC.SCoV.over.Mock", "adj.P.Val.SCoV.over.Mock")]
colnames(top1schmidt) <- c("gene", "Schmidt_FC", "Schmidt_pvalue")


# Import and process Klann
all_klann <- readxl::read_xlsx("../data/Klann_TS2.xlsx", sheet = 1) %>% data.frame() %>% 
  dplyr::select(c("Gene.Symbol", "log.2.Ratio.Virus.Mock","adjusted.P.value"))  %>% arrange((adjusted.P.value))

top1klann <- all_klann[!duplicated(all_klann$Gene.Symbol),]
colnames(top1klann) <- c("gene", "klann_FC", "klann_pvalue")

# Merge altogether
mdf <- merge(top1schmidt, top1klann, by = "gene")

cor(mdf[mdf$Schmidt_pvalue < 0.01 & mdf$klann_pvalue < 0.01 ,c(2,4)], use = "pairwise.complete.obs")
ss <- mdf %>% filter(Schmidt_pvalue < 0.01 & klann_pvalue < 0.01)
table(ss$Schmidt_FC > 0, ss$klann_FC > 0)

# Merge sipmle
p1 <- ggplot(mdf %>% filter(Schmidt_pvalue < 0.01 & klann_pvalue < 0.01),
       aes(x = Schmidt_FC, y = klann_FC)) +
  geom_point(size = 0.3) + pretty_plot(fontsize = 8) +
  labs(x = "log2 FC Schmidt", y = "log2 FC Klann") +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_hline(yintercept = 0, linetype = 2)
cowplot::ggsave2(p1, file = "../output/Klann_bulkProteome_viz.pdf", 
                 width = 2, height = 2)
#---

# Import and process bojkova
all_bojkova <- readxl::read_xlsx("../data/Bojkova_32408336_TableS2.xlsx") %>% data.frame() %>% 
  dplyr::select(c("Gene.Symbol", "Ratio.24h","P.value.24h"))  %>% arrange((P.value.24h))
top1boj <- all_bojkova[!duplicated(all_bojkova$Gene.Symbol),]
colnames(top1boj) <- c("gene", "Bojkova_FC", "Bojkova_pvalue")

# Import and process Bouhaddou
all_bouhaddou <- readxl::read_xlsx("../data/Bouhaddou2020_TS1.xlsx") %>% data.frame() %>% 
  dplyr::select(c("Gene_Name", "Inf_24Hr.log2FC","Inf_24Hr.adj.pvalue"))  %>% arrange((Inf_24Hr.adj.pvalue))
top1bouhaddou <- all_bouhaddou[!duplicated(all_bouhaddou$Gene_Name),]
colnames(top1bouhaddou) <- c("gene", "Bouhaddou_FC", "Bouhaddou_pvalue")
top1bouhaddou$Bouhaddou_FC <- as.numeric(top1bouhaddou$Bouhaddou_FC )

