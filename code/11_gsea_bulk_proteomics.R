library(fgsea)
library(dplyr)
library(data.table)
library(org.Hs.eg.db)
library(clusterProfiler)
library(msigdbr)


viral_genes <- c("N", "S", "3a", "M", "ORF9b", "PP1AB")

# Import differential gene expression
all <- fread("../data/all_proteomics.txt") %>% arrange(desc(Log.P.Value.SCoV.over.Mock)) %>%
  filter(!(geneSymbol %in% viral_genes))
diff_df <- all[!duplicated(all$geneSymbol),]

tl <- bitr(diff_df$gene, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db", drop = FALSE) 
tl <- tl[!duplicated(tl),]; name_vec <- tl[[2]]; names(name_vec) <- tl[[1]]
diff_df$gene_id <- name_vec[as.character(diff_df$gene)]
diff_df <- diff_df[complete.cases(diff_df),]

# Do GSEA
geneList <- diff_df[["logFC.SCoV.over.Mock"]] * diff_df[["Log.P.Value.SCoV.over.Mock"]]; names(geneList) <- diff_df[["geneSymbol"]]
geneList <- sort(geneList, decreasing = FALSE)

h_df = msigdbr(species = "Homo sapiens") %>% dplyr::filter(gs_cat == "C5" & gs_subcat == "BP")
pathway_list = h_df %>% split(x = .$gene_symbol, f = .$gs_name)

fgseaRes <- fgseaMultilevel(pathways = pathway_list, 
                            stats    = geneList,
                            minSize  = 15,
                            maxSize  = 500)

if(FALSE){
  collapsedPathways <- collapsePathways(fgseaRes[order(pval)][padj < 0.01], 
                                        pathway_list, geneList)
  mainPathways <- fgseaRes[pathway %in% collapsedPathways$mainPathways][
    order(-abs(NES)), pathway]
  plotGseaTable(pathway_list[mainPathways], geneList, fgseaRes, 
                gseaParam = 0.5)
}
fgseaRes[,c(1:7)] %>% arrange((padj)) %>% head(100)


fgseaRes[,c(1:7)] %>% arrange((padj)) %>% filter(pathway %in% c("GO_RESPONSE_TO_TYPE_I_INTERFERON",
                                                                "GO_REGULATION_OF_MAPK_CASCADE", "GO_POSITIVE_REGULATION_OF_MAP_KINASE_ACTIVITY"))


