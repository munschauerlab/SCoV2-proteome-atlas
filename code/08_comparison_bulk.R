library(data.table)
library(dplyr)
library(Matrix)
library(BuenColors)

# Import gene names
all <- fread("../data/SCoV2_bulkProteome.txt") %>% arrange(log(P.Value.SCoV.over.Mock)) %>%
  filter(species == "HOMO SAPIENS") 

top1ms <- all[!duplicated(all$geneSymbol) & all$adj.P.Val.SCoV.over.Mock < 0.1,]
top1schmidt <- top1ms[,c("geneSymbol", "logFC.SCoV.over.Mock", "adj.P.Val.SCoV.over.Mock")]
colnames(top1schmidt) <- c("gene", "Schmidt_FC", "Schmidt_pvalue")

# Import and process bojkova
all_bojkova <- readxl::read_xlsx("../data/Bojkova_32408336_TableS2.xlsx") %>% data.frame() %>% 
  dplyr::select(c("Gene.Symbol", "Ratio.24h","P.value.24h"))  %>% arrange((P.value.24h))
top1boj <- all_bojkova[!duplicated(all_bojkova$Gene.Symbol),]
colnames(top1boj) <- c("gene", "Bojkova_FC", "Bojkova_pvalue")

# Import and process bojkova
all_bouhaddou <- readxl::read_xlsx("../data/Bouhaddou2020_TS1.xlsx") %>% data.frame() %>% 
  dplyr::select(c("Gene_Name", "Inf_24Hr.log2FC","Inf_24Hr.adj.pvalue"))  %>% arrange((Inf_24Hr.adj.pvalue))
top1bouhaddou <- all_bouhaddou[!duplicated(all_bouhaddou$Gene_Name),]
colnames(top1bouhaddou) <- c("gene", "Bouhaddou_FC", "Bouhaddou_pvalue")
top1bouhaddou$Bouhaddou_FC <- as.numeric(top1bouhaddou$Bouhaddou_FC )


# Merge altogether
mdf <- merge(merge(top1schmidt, top1boj, by = "gene"),top1bouhaddou, by = "gene")

cor(mdf[,c(2,4,6)], use = "pairwise.complete.obs")

# Merge sipmle
qplot(merge(top1schmidt, top1boj, by = "gene")[[2]],
      merge(top1schmidt, top1boj, by = "gene")[[4]]) + 
  labs(x = "Schmidt FC", y = "Bojkova FC")
