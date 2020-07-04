library(fgsea)
library(dplyr)
library(data.table)
library(org.Hs.eg.db)
library(clusterProfiler)
library(msigdbr)
library(xlsx)

"%ni%" <- Negate("%in%")

# Import differential gene expression
all <- fread("../data/SCoV2_bulkProteome.txt") %>% arrange(log(P.Value.SCoV.over.Mock)) %>%
  filter(species == "HOMO SAPIENS") 
diff_df <- all[!duplicated(all$geneSymbol),]

tl <- bitr(diff_df$gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db", drop = FALSE) 
tl <- tl[!duplicated(tl),]; name_vec <- tl[[2]]; names(name_vec) <- tl[[1]]
diff_df$gene_id <- name_vec[as.character(diff_df$gene)]
diff_df <- diff_df[complete.cases(diff_df),]

# Do GSEA
geneList <- diff_df[["logFC.SCoV.over.Mock"]] * -1*log10(diff_df[["P.Value.SCoV.over.Mock"]]); names(geneList) <- diff_df[["geneSymbol"]]
geneList <- sort(geneList, decreasing = FALSE)

h_df_hallmark = msigdbr(species = "Homo sapiens") %>% dplyr::filter(gs_cat == "H")
pathway_list_hallmark = h_df_hallmark %>% split(x = .$gene_symbol, f = .$gs_name)
names(pathway_list_hallmark) <- gsub("HALLMARK_", "", names(pathway_list_hallmark))

h_df_c5bp = msigdbr(species = "Homo sapiens") %>% dplyr::filter(gs_cat == "C5" & gs_subcat == "BP")
pathway_list_c5bp = h_df_c5bp %>% split(x = .$gene_symbol, f = .$gs_name)


fgseaRes_hallmark <- fgseaMultilevel(pathways = pathway_list_hallmark, 
                                     stats    = geneList,
                                     minSize  = 15,
                                     maxSize  = 500)

fgseaRes_hallmark[,c(1:7)] %>% arrange(pval)

fgseaRes_c5bp <- fgseaMultilevel(pathways = pathway_list_c5bp, 
                                 stats    = geneList,
                                 minSize  = 15,
                                 maxSize  = 500)


# Visualize the hallmark pathways in figure
collapsedPathways <- collapsePathways(fgseaRes_hallmark[order(pval)][padj < 0.1], 
                                      pathway_list_hallmark, geneList)
mainPathways <- fgseaRes_hallmark[pathway %in% names(collapsedPathways$parentPathways)][
  order(-abs(NES)*(sign(NES) + 2)), pathway]

selected <- c("TNFA_SIGNALING_VIA_NFKB", "TGF_BETA_SIGNALING", "IL6_JAK_STAT3_SIGNALING", "INTERFERON_GAMMA_RESPONSE", "INTERFERON_ALPHA_RESPONSE")
others <- mainPathways[mainPathways %ni% selected]

pdf("../output/GSEA_hallmark_viz-main.pdf", width = 6, height = 3)
plotGseaTable(pathway_list_hallmark[selected], geneList, fgseaRes_hallmark, 
              gseaParam = 0.5)
dev.off()


pdf("../output/GSEA_hallmark_viz-supp.pdf", width = 6, height = 5)
plotGseaTable(pathway_list_hallmark[others], geneList, fgseaRes_hallmark, 
              gseaParam = 0.5)
dev.off()

# Export pathways
selected_pathways <- c("GO_RESPONSE_TO_TYPE_I_INTERFERON","GO_REGULATION_OF_MAPK_CASCADE", "GO_POSITIVE_REGULATION_OF_MAP_KINASE_ACTIVITY")
fgseaRes_c5bp %>% filter(pathway %in% selected_pathways)

fgseaRes_hallmark[,c(1:7)] %>% arrange((pval)) %>% 
  write.xlsx("../output/TableS6_GO_BP_Proteome_GSEA.xlsx", sheetName = "Hallmark", 
             col.names = TRUE, row.names = FALSE, append = TRUE)
fgseaRes_c5bp[,c(1:7)] %>% arrange((padj)) %>% 
  write.xlsx("../output/TableS6_GO_BP_Proteome_GSEA.xlsx", sheetName = "MSigDB C5 BP", 
             col.names = TRUE, row.names = FALSE, append = TRUE)




