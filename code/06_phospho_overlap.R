library(data.table)
library(readxl)
library(dplyr)

Bouhaddou <- readxl::read_xlsx("../data/Bouhaddou2020_TS1.xlsx") %>%
  data.frame() %>%
  filter(Inf_00Hr.adj.pvalue < 0.05 | Inf_02Hr.adj.pvalue < 0.05 | Inf_04Hr.adj.pvalue < 0.05 | Inf_08Hr.adj.pvalue < 0.05 | Inf_12Hr.adj.pvalue < 0.05 |  Inf_24Hr.adj.pvalue < 0.05 )

Klann <- readxl::read_xlsx("../data/Klann_TS1.xlsx", sheet = 2) %>%
  data.frame() %>%
  filter(adjusted.P < 0.5)

Stukalov <- fread("../data/Stukalov_TableS7.txt") %>%
  data.frame() %>%
  filter(q.value.1 < 0.05 | q.value < 0.05)

rap <- fread("../data/SCoV2_RAPms.txt") %>%
  filter(adj.P.Val.SCoV2.over.RMRP < 0.2 & species == "HOMO SAPIENS" & logFC.SCoV2.over.RMRP > 0)

table(rap$geneSymbol %in% Bouhaddou$Gene_Name)
table(rap$geneSymbol %in% Klann$T..Gene.name)
table(rap$geneSymbol %in% Stukalov$gene.name)

intersect(rap$geneSymbol[rap$geneSymbol %in% Bouhaddou$Gene_Name],
rap$geneSymbol[rap$geneSymbol %in% Stukalov$gene.name])
