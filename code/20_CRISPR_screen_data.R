library(readxl)
library(dplyr)
library(data.table)

gene_zscores <- readxl::read_xlsx("../data/Wei_Cell2020_TableS1.xlsx") %>% data.frame()
label <- colnames(gene_zscores)
v2_avg <- gene_zscores[,c(1,5,6,7,8,9)] %>% 
  reshape2::melt(id.vars = c("Gene")) %>%
  mutate(p_value = 2*pnorm(-abs(value))) %>%
  group_by(Gene) %>%
  summarise(mean_z = mean(value),
            combined_pvalue = pchisq(-2 * sum(log(p_value)),df=2*length(p_value),lower=FALSE)) %>% # Fisher's method
  mutate(fdr = p.adjust(combined_pvalue, method = 'BH')) 

# Annotate with RAP
rap <- fread("../data/SCoV2_RAPms.txt") %>%
  filter(adj.P.Val.SCoV2.over.RMRP < 0.2 & species == "HOMO SAPIENS" & logFC.SCoV2.over.RMRP > 0)
rap_genes <- rap %>% pull(geneSymbol)

df <- data.frame(v2_avg)
df$in_RAP <- df$Gene %in% rap_genes
table(df$in_RAP)

rap_genes[!(rap_genes %in% df$Gene)]
xlsx::write.xlsx(df, "../output/Wilen_CRISPRscreen_FDR_RAP.xlsx", sheetName = "Sheet1",
            col.names = TRUE, row.names = FALSE, append = FALSE)
